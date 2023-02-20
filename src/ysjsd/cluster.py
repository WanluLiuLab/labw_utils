from __future__ import annotations

import json
import logging
import multiprocessing
import os
import signal
import subprocess
import threading
import time
from typing import Union, Optional, Dict, List

import psutil
import sqlalchemy
from sqlalchemy.orm import sessionmaker

from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper import logger_helper
from libysjs.cluster import YSJSDLoad
from libysjs.job import CoreYSJSJob
from libysjs.submission import YSJSSubmission
from ysjsd.orm import SQLAlchemyDeclarativeBase, ServerSideYSJSDConfigTable, YSJSDVersionTable, SubmissionTable
from ysjsd.config import ServerSideYSJSDConfig


class YSJSD(threading.Thread):
    __submission_queue_lock: threading.Lock
    _db_write_lock: threading.Lock
    _config: ServerSideYSJSDConfig
    _current_cpu: Union[int, float]
    _current_mem: Union[int, float]
    _should_terminate: bool
    _schedule_method: str
    _lh: logging.Logger
    _dbe: sqlalchemy.engine.Engine
    _submissions_waiting_to_be_proceed: List[YSJSSubmission]

    def __init__(self, config: ServerSideYSJSDConfig):
        super().__init__()
        self.__submission_queue_lock = threading.Lock()
        self._db_write_lock = threading.Lock()
        self._config = config
        os.makedirs(self._config.var_directory, exist_ok=True)
        os.makedirs(os.path.join(self._config.var_directory, "submission"), exist_ok=True)

        log_file_path = os.path.join(self._config.var_directory, "ysjsd.log")
        self._lh = logger_helper.get_logger(
            name="YSJSD",
            level=logger_helper.TRACE,
            log_file_name=log_file_path,
            log_file_level=logger_helper.TRACE
        )
        self._lh.info("Logger set up at %s", log_file_path)
        self._config.validate()
        self._should_terminate = False
        self._current_cpu = self._config.total_cpu
        self._current_mem = self._config.total_mem
        self._schedule_method = self._config.schedule_method
        self._submissions_waiting_to_be_proceed = []
        dburl = f"sqlite:///{self._config.var_directory}/foo.db"
        with self._db_write_lock:
            self._dbe = sqlalchemy.engine.create_engine(
                url=dburl
            )
            metadata = SQLAlchemyDeclarativeBase.metadata
            create_drop_params = {"bind": self._dbe, "checkfirst":True}
            metadata.tables[ServerSideYSJSDConfigTable.__tablename__].drop(**create_drop_params)
            metadata.tables[YSJSDVersionTable.__tablename__].drop(**create_drop_params)
            metadata.create_all(**create_drop_params)
            with sessionmaker(bind=self._dbe)() as session:
                session.add(
                    ServerSideYSJSDConfigTable(**self._config.to_dict())
                )
                for name, version in ServerSideYSJSDConfig.dump_versions().items():
                    session.add(YSJSDVersionTable(name=name, version=version))
                session.commit()
        self._lh.info("Initialized database %s", dburl)
        self._lh.info("Initialization Finished")

    def receive_submission(self, submission: YSJSSubmission):
        sid = submission.submission_id
        with get_writer(
                os.path.join(self._config.var_directory, "submission", f"{sid}.json")
        ) as writer:
            json.dump(submission.to_dict(), writer)
        with self._db_write_lock:
            with sessionmaker(bind=self._dbe)() as session:
                session.add(
                    SubmissionTable(**submission.to_dict())
                )
                session.commit()
        with self.__submission_queue_lock:
            self._submissions_waiting_to_be_proceed.append(submission)
    
    def get_submission_aggressive(self, avail_cpu: float, avail_mem: float) -> Optional[YSJSSubmission]:
        if not self._submissions_waiting_to_be_proceed:
            return None
        with self.__submission_queue_lock:
            for i in range(len(self._submissions_waiting_to_be_proceed)):
                submission = self._submissions_waiting_to_be_proceed[i]
                if submission.cpu < avail_cpu and submission.mem < avail_mem:
                    return self._submissions_waiting_to_be_proceed.pop(i)
        return None


    def get_submission_fifo(self, avail_cpu: float, avail_mem: float) -> Optional[YSJSSubmission]:
        if not self._submissions_waiting_to_be_proceed:
            return None
        with self.__submission_queue_lock:
            submission = self._submissions_waiting_to_be_proceed[0]
            if submission.cpu < avail_cpu and submission.mem < avail_mem:
                return self._submissions_waiting_to_be_proceed.pop(0)
        return None
    

    def run(self):
        self._lh.info("Started at http://localhost:%s", self._config.ysjs_port)

        while not self._should_terminate:
            time.sleep(1)
        self._lh.info("Terminated")

    def terminate(self):
        self._lh.info("Received termination signal")
        self._should_terminate = True

    @property
    def real_load(self) -> YSJSDLoad:
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
