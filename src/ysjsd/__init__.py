from __future__ import annotations

import json
import logging
import multiprocessing
import os
import platform
import signal
import subprocess
import sys
import threading
import time
from typing import Union, Dict, Mapping, Any, Optional, List

import flask
import psutil

import labw_utils
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.serializer.toml import AbstractTOMLSerializable
from labw_utils.commonutils.stdlib_helper import logger_helper
from libysjs import __version__
from libysjs.cluster import YSJSDLoad, CoreYSJSDConfig, AVAILABLE_SCHEDULING_METHOD
from libysjs.job import CoreYSJSJob
from libysjs.queue import YSJSDQueueConfig
from libysjs.submission import YSJSSubmission


class YSJSD(threading.Thread):
    _lock: threading.Lock
    _config: ServerSideYSJSDConfig
    _current_cpu: Union[int, float]
    _current_mem: Union[int, float]
    _should_terminate: bool
    _schedule_method: str
    _lh: logging.Logger
    _queue: Dict[str, ServerSideYSJSQueue]

    def __init__(self, config: ServerSideYSJSDConfig):
        super().__init__()
        self._lock = threading.Lock()
        self._config = config
        os.makedirs(self._config.var_directory)
        os.makedirs(os.path.join(self._config.var_directory, "submission"))

        self._lh = logger_helper.get_logger(
            name="YSJSD",
            level=logger_helper.TRACE,
            log_file_name=os.path.join(self._config.var_directory, "ysjsd.log"),
            log_file_level=logger_helper.TRACE
        )
        self._lh.info("Logger set up")
        self._config.validate()
        self._should_terminate = False
        self._current_cpu = self._config.total_cpu
        self._current_mem = self._config.total_mem
        self._schedule_method = self._config.schedule_method
        self._lh.info("Initialization Finished")

    def receive_submission(self, submission: YSJSSubmission):
        sid = submission.submission_id
        with get_writer(
                os.path.join(self._config.var_directory, "submission", f"{sid}.json")
        ) as writer:
            json.dump(submission.to_dict(), writer)

    def run(self):
        self._lh.info("Started at http://localhost:%s", self._config.ysjs_port)
        while not self._should_terminate:
            time.sleep(1)
        self._lh.info("Terminated")

    def terminate(self):
        self._lh.info("Received termination signal")
        self._should_terminate = True

    @staticmethod
    def get_real_load() -> YSJSDLoad:
        real_total_cpu = multiprocessing.cpu_count()
        return YSJSDLoad(
            real_total_cpu=real_total_cpu,
            real_avail_cpu=(1 - psutil.cpu_percent(1)) * real_total_cpu,
            real_total_mem=psutil.virtual_memory().total,
            real_avail_mem=psutil.virtual_memory().available
        )


class ServerSideYSJSJob(CoreYSJSJob, threading.Thread):
    _p: Optional[subprocess.Popen]
    _status: str

    def __init__(self, submission: YSJSSubmission):
        super().__init__()
        self._submission = submission
        self._p = None
        self._status = "pending"
        self._retv = None

    def run(self):
        if self._status == "canceled":
            return
        self._status = "starting"
        if self._submission.stdin is None:
            stdin = subprocess.DEVNULL
        else:
            stdin = open(self._submission.stdin, "rb")
        if self._submission.stdout is None:
            stdout = subprocess.DEVNULL
        else:
            stdout = open(self._submission.stdout, "wb")
        if self._submission.stderr is None:
            stderr = subprocess.DEVNULL
        else:
            stderr = open(self._submission.stdin, "wb")

        self._p = subprocess.Popen(
            args=[
                self._submission.shell_path,
                self._submission.script_path
            ],
            cwd=self._submission.cwd,
            env=self._submission.env,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            close_fds=True
        )

        self._status = "started"
        self._retv = self._p.wait()
        self._status = "finished" if self._retv == 0 else "failed"

    def send_signal(self, _signal: int):
        self._p.send_signal(_signal)

    @property
    def pid(self) -> Optional[int]:
        if self._p is not None:
            return self._p.pid
        else:
            return None

    def cancel(self):
        self._status = "canceled"

    @property
    def status(self) -> str:
        return self._status

    def kill(self, timeout: float):
        """
        Recursively terminate a process tree
        """
        try:
            p = psutil.Process(self._p.pid)
            children = p.children(recursive=True)
        except psutil.Error:
            return
        processes_to_be_killed = [p, *children]
        for process_to_be_killed in processes_to_be_killed:
            try:
                process_to_be_killed.send_signal(signal.SIGTERM)
            except psutil.Error:
                pass
        time.sleep(timeout)
        try:
            p = psutil.Process(self._p.pid)
            children = p.children(recursive=True)
        except psutil.Error:
            return
        processes_to_be_killed = [p, *children]
        for process_to_be_killed in processes_to_be_killed:
            try:
                process_to_be_killed.send_signal(signal.SIGKILL)
            except psutil.Error:
                pass


class ServerSideWorkbench(threading.Thread):
    _max_concurrent_jobs: int
    _jobs: Dict[str, ServerSideYSJSJob]
    _should_terminate: bool
    _kill_timeout: float

    def terminate(self):
        self._should_terminate = True

    def __init__(
            self,
            max_concurrent_jobs: int,
            kill_timeout: float
    ):
        super().__init__()
        self._should_terminate = False
        self._max_concurrent_jobs = max_concurrent_jobs
        self._kill_timeout = kill_timeout

    def run(self):
        ...

    def job_send_signal(self, submission_id: str, _signal: int):
        self._jobs[submission_id].send_signal(_signal)

    def job_kill(self, submission_id: str):
        self._jobs[submission_id].kill(self._kill_timeout)

    def queue_kill(self, queue_name: str):
        kill_thread_pool = []
        for job in self._jobs.values():
            if job.submission.queue_name == queue_name:
                kill_thread_pool.append(
                    threading.Thread(target=job.kill, args=[self._kill_timeout])
                )
                kill_thread_pool[-1].start()
        for kill_thread in kill_thread_pool:
            kill_thread.join()

    def queue_send_signal(self, queue_name: str, _signal: int):
        for job in self._jobs.values():
            if job.submission.queue_name == queue_name:
                job.send_signal(_signal)


class ServerSideYSJSDConfig(CoreYSJSDConfig, AbstractTOMLSerializable):
    @staticmethod
    def _validate_versions(versions: Mapping[str, Any]) -> None:
        _lh = logger_helper.get_logger("YSJSD - LOADING CONFIG")
        compile_time_version_dict = dict(versions)
        run_time_version_dict = ServerSideYSJSDConfig._dump_versions()
        for version_key, run_time_version_value in run_time_version_dict.items():
            try:
                compile_time_version_value = compile_time_version_dict.pop(version_key)
            except KeyError:
                compile_time_version_value = None
            if compile_time_version_value != run_time_version_value:
                _lh.warning(
                    "Package %s have different version information: Compile (%s) != Run (%s)",
                    version_key, compile_time_version_value, run_time_version_value
                )
        for remaining_compile_time_version_key, remaining_compile_time_version_value in compile_time_version_dict.items():
            _lh.warning(
                "Package %s have different version information: Compile (%s) != Run (%s)",
                remaining_compile_time_version_key, remaining_compile_time_version_value, None
            )

    @staticmethod
    def _dump_versions() -> Mapping[str, Any]:
        return {
            "labw_utils": labw_utils.__version__,
            "ysjsd": __version__,
            "flask": flask.__version__,
            "python": ".".join(map(str, (sys.version_info.major, sys.version_info.minor, sys.version_info.micro)))
        }

    @staticmethod
    def _dump_metadata() -> Mapping[str, Any]:
        return {
            "time_gmt": time.asctime(time.gmtime()),
            "time_local": time.asctime(time.localtime()),
            "platform": platform.uname()._asdict()
        }

    @classmethod
    def new(cls, config_file_path: str):
        return cls(
            name="ylsjs_cluster",
            description="No description",
            ysjs_port="8080",
            var_directory=os.path.join(os.path.abspath("."), "var"),
            config_file_path=os.path.abspath(config_file_path),
            total_cpu=multiprocessing.cpu_count(),
            total_mem=(psutil.virtual_memory().total + psutil.swap_memory().total) * 0.8,
            schedule_method="FIFO",
            queue_conf=[
                YSJSDQueueConfig("root", "The root queue", 1024)
            ]
        )

    def validate(self):
        _lh = logger_helper.get_logger("YSJSD")
        if self.schedule_method not in AVAILABLE_SCHEDULING_METHOD:
            raise ValueError(
                f"Illegal scheule_method {self.schedule_method}, "
                f"should be in {str(AVAILABLE_SCHEDULING_METHOD)}"
            )
        max_total_cpu = multiprocessing.cpu_count()
        if self.total_cpu > max_total_cpu:
            _lh.warning(
                "Configured CPU number %.2f larger than total CPU number %d",
                self.total_cpu,
                max_total_cpu
            )
        recommended_total_mem = (psutil.virtual_memory().total + psutil.swap_memory().total) * 0.8
        if self.total_mem > recommended_total_mem:
            _lh.warning(
                "Configured memory size %.2f larger than total CPU number %.2f",
                self.total_mem,
                recommended_total_mem
            )


class ServerSideYSJSQueue:
    _submissions: List[YSJSSubmission]
    _lock: threading.Lock
    _queue_config: YSJSDQueueConfig

    def __init__(
            self,
            queue_config: YSJSDQueueConfig
    ):
        self._lock = threading.Lock()
        self._queue_config = queue_config
        self._submissions = []

    def add(self, submission: YSJSSubmission):
        with self._lock:
            self._submissions.append(submission)

    def get_submission_aggressive(
            self,
            avail_cpu: Union[float, int],
            avail_mem: Union[float, int]
    ) -> Optional[YSJSSubmission]:
        with self._lock:
            for submission in self._submissions:
                if submission.cpu < avail_cpu and submission.mem < avail_mem:
                    return submission
            return None

    def get_submission_fifo(
            self,
            avail_cpu: Union[float, int],
            avail_mem: Union[float, int]
    ) -> Optional[YSJSSubmission]:
        with self._lock:
            if not self._submissions:
                return None
            else:
                submission = self._submissions[0]
                if submission.cpu < avail_cpu and submission.mem < avail_mem:
                    self._submissions.pop(0)
                    return submission

    @property
    def earliest_submission_time(self) -> Optional[float]:
        with self._lock:
            if not self._submissions:
                return None
            else:
                return self._submissions[0].submission_time

    @property
    def __len__(self):
        return len(self._submissions)
