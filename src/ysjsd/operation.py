from __future__ import annotations

import json
import logging
import multiprocessing
import os
import sys
import threading
import time
from typing import Union, Optional, Dict, Iterable

import psutil
import sqlalchemy
from sqlalchemy.orm import sessionmaker

from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.shell_utils import rm_rf
from labw_utils.commonutils.stdlib_helper import logger_helper
from libysjs.ds.ysjs_submission import YSJSSubmission
from libysjs.operation import YSJSDLoad
from ysjsd.ds.ysjs_job import ServerSideYSJSJob
from ysjsd.ds.ysjsd_config import ServerSideYSJSDConfig
from ysjsd.orm import SQLAlchemyDeclarativeBase
from ysjsd.orm.ysjs_submission_table import YSJSSubmissionTable
from ysjsd.orm.ysjsd_config_table import ServerSideYSJSDConfigTable
from ysjsd.orm.ysjsd_version_table import YSJSDVersionTable


class YSJSDException(RuntimeError):
    ...


class JobNotExistException(YSJSDException):
    ...


class IllegalOperationException(YSJSDException):
    ...


class NotAvailableException(YSJSDException):
    ...


def get_job(job_id: int) -> ServerSideYSJSJob:
    return ...


class YSJSD(threading.Thread):
    _job_queue_lock: threading.Lock
    _db_write_lock: threading.Lock
    _config: ServerSideYSJSDConfig
    _current_cpu: Union[int, float]
    _current_mem: Union[int, float]
    _state: str
    _schedule_method: str
    _lh: logging.Logger
    _dbe: sqlalchemy.engine.Engine
    _job_queue_pending: Dict[int, ServerSideYSJSJob]
    _job_queue_running: Dict[int, ServerSideYSJSJob]
    _latest_job_id: int
    _last_job_id_filename: str
    _lock_filename: str

    def __init__(self, config: ServerSideYSJSDConfig):
        super().__init__()
        self._job_queue_lock = threading.Lock()
        self._db_write_lock = threading.Lock()
        self._config = config

        # Create log
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

        # Process PID lock
        self._lock_filename = os.path.join(self._config.var_directory, "pid.lock")
        try:
            with get_reader(self._lock_filename) as last_job_id_reader:
                prev_pid = int(last_job_id_reader.read())
                if psutil.pid_exists(prev_pid):
                    self._lh.error(
                        "Previous lock %s (pid=%d) still running!",
                        self._lock_filename,
                        prev_pid
                    )
                    sys.exit(1)
                else:
                    self._lh.warning(
                        "Previous lock %s (pid=%d) invalid; will be removed",
                        self._lock_filename,
                        prev_pid
                    )
                    rm_rf(self._lock_filename)
        except (ValueError, FileNotFoundError):
            self._lh.warning("Previous lock %s invalid; will be removed", self._lock_filename)
            rm_rf(self._lock_filename)
        with get_writer(self._lock_filename) as lock_writer:
            lock_writer.write(f"{os.getpid()}\n")

        # Other configs
        self._config.validate()
        self._state = "starting"
        self._current_cpu = self._config.total_cpu
        self._current_mem = self._config.total_mem
        self._schedule_method = self._config.schedule_method
        self._job_queue_pending = {}
        self._job_queue_running = {}

        # Connect to DB
        dburl = f"sqlite:///{self._config.var_directory}/foo.db"
        with self._db_write_lock:
            self._dbe = sqlalchemy.engine.create_engine(
                url=dburl
            )
            metadata = SQLAlchemyDeclarativeBase.metadata
            create_drop_params = {"bind": self._dbe, "checkfirst": True}
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

        # Load last job_id
        self._latest_job_id = 0
        self._last_job_id_filename = os.path.join(self._config.var_directory, "last_job.txt")
        try:
            with get_reader(self._last_job_id_filename) as last_job_id_reader:
                self._latest_job_id = int(last_job_id_reader.read())
        except (ValueError, FileNotFoundError):
            self._lh.warning("Previous last_job_id file %s invalid; will be removed", self._last_job_id_filename)
            rm_rf(self._last_job_id_filename)
        # TODO: Load from Database

        # Finished
        self._lh.info("Initialization Finished")

    def receive_submission(self, submission: YSJSSubmission) -> int:
        """
        Read one submission and load it into pending job queue,
        and return its jobid
        """
        if self._state != "running":
            raise NotAvailableException
        sid = submission.submission_id
        with get_writer(
                os.path.join(self._config.var_directory, "submission", f"{sid}.json")
        ) as writer:
            json.dump(submission.to_dict(), writer)
        with self._db_write_lock:
            with sessionmaker(bind=self._dbe)() as session:
                session.add(
                    YSJSSubmissionTable(**submission.to_dict())
                )
                session.commit()
        with self._job_queue_lock:
            new_job = ServerSideYSJSJob.new(
                dbe= self._dbe,
                db_write_lock=self._db_write_lock,
                submission=submission,
                job_id=self._latest_job_id
            )
            self._job_queue_pending[self._latest_job_id] = new_job
            reti = self._latest_job_id
            self._latest_job_id += 1
            with get_writer(self._last_job_id_filename) as last_job_id_writer:
                last_job_id_writer.write(f"{self._latest_job_id}\n")
        return reti

    def _fetch_pending_job(self, avail_cpu: float, avail_mem: float) -> Optional[ServerSideYSJSJob]:
        if not self._job_queue_pending or self._state != "running":
            return None
        with self._job_queue_lock:
            if self._config.schedule_method == "FIFO":
                try:
                    job = self._job_queue_pending[0]
                    if job.submission.cpu < avail_cpu and job.submission.mem < avail_mem:
                        return self._job_queue_pending.pop(0)
                except KeyError:
                    return None
            elif self._config.schedule_method == "AGGRESSIVE":
                for i in range(len(self._job_queue_pending)):
                    job = self._job_queue_pending[i]
                    if job.submission.cpu < avail_cpu and job.submission.mem < avail_mem:
                        return self._job_queue_pending.pop(i)
        return None

    def job_cancel(self, job_id: int):
        if self._state == "starting":
            raise NotAvailableException
        with self._job_queue_lock:
            try:
                job_to_cancel = self._job_queue_pending.pop(job_id)
            except KeyError as e:
                raise JobNotExistException from e
        job_to_cancel.cancel()

    def job_send_signal(self, job_id: int, _signal: int):
        if self._state == "starting":
            raise NotAvailableException
        try:
            self._job_queue_running[job_id].send_signal(_signal)
        except KeyError as e:
            raise JobNotExistException from e

    def job_kill(self, job_id: int):
        if self._state == "starting":
            raise NotAvailableException
        try:
            self._job_queue_running[job_id].kill(self._config.kill_timeout)
        except KeyError as e:
            raise JobNotExistException from e

    def query(self, ) -> Iterable[int]:
        ...

    def apply(
            self,
            job_ids: Iterable[int],
            operation: str,
            **extra_params
    ):
        if self._state == "starting":
            raise NotAvailableException
        if operation == "cancel":
            map(self.job_cancel, job_ids)
        elif operation == "kill":
            map(self.job_kill, job_ids)
        elif operation == "send_signal":
            map(lambda job_id: self.job_send_signal(job_id=job_id, **extra_params), job_ids)
        else:
            raise IllegalOperationException

    def run(self):
        self._lh.info("Started at http://localhost:%s", self._config.ysjs_port)
        self._state = "running"
        while self._state != "terminating":
            if len(self._job_queue_running) < self._config.max_concurrent_jobs:
                job_fetched = self._fetch_pending_job(avail_cpu=self._current_cpu, avail_mem=self._current_mem)
                if job_fetched is not None:
                    self._current_cpu -= job_fetched.submission.cpu
                    self._current_mem -= job_fetched.submission.mem
                    self._job_queue_running[job_fetched.job_id] = job_fetched
                    job_fetched.start()
            for job_id in list(self._job_queue_running.keys()):
                job = self._job_queue_running[job_id]
                if job.status == "finished":
                    self._job_queue_running.pop(job_id)
                    self._current_cpu += job.submission.cpu
                    self._current_mem += job.submission.mem
            time.sleep(0.1)

        self._lh.info("Terminating")
        self.apply(
            self._job_queue_pending.keys(),
            operation="cancel"
        )
        self.apply(
            self._job_queue_running.keys(),
            operation="kill"
        )
        rm_rf(self._lock_filename)
        self._lh.info("Terminated")

    def terminate(self):
        self._lh.info("Received termination signal")
        self._state = "terminating"

    @property
    def real_load(self) -> YSJSDLoad:
        real_total_cpu = multiprocessing.cpu_count()
        return YSJSDLoad(
            real_total_cpu=real_total_cpu,
            real_avail_cpu=(1 - psutil.cpu_percent(1) / 100) * real_total_cpu,
            real_total_mem=psutil.virtual_memory().total,
            real_avail_mem=psutil.virtual_memory().available
        )