from __future__ import annotations

import json
import logging
import multiprocessing
import os
import sys
import threading
import time

import psutil
import sqlalchemy
from sqlalchemy.orm import sessionmaker

from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.stdlib_helper import logger_helper
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
from labw_utils.typing_importer import Union, Optional, Dict, Iterable
from libysjs.ds.ysjs_submission import YSJSSubmission
from libysjs.ds.ysjsd_status import YSJSDStatus
from libysjs.operation import YSJSDLoad
from ysjsd.ds.ysjs_job import ServerSideYSJSJob
from ysjsd.ds.ysjsd_config import ServerSideYSJSDConfig
from ysjsd.orm import SQLAlchemyDeclarativeBase
from ysjsd.orm.ysjs_job_table import YSJSJobTable
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
    _start_time: float

    def __init__(self, config: ServerSideYSJSDConfig):
        super().__init__()
        self._start_time = 0
        self._job_queue_lock = threading.Lock()
        self._db_write_lock = threading.Lock()
        self._config = config

        # Create log
        os.makedirs(self._config.var_directory_path, exist_ok=True)
        os.makedirs(os.path.join(self._config.var_directory_path, "submission"), exist_ok=True)
        log_file_path = os.path.join(self._config.var_directory_path, "ysjsd.log")
        self._lh = logger_helper.get_logger(
            name="YSJSD",
            level=logger_helper.TRACE,
            log_file_name=log_file_path,
            log_file_level=logger_helper.TRACE
        )
        self._lh.info("Logger set up at %s", log_file_path)

        # Process PID lock
        self._lock_filename = os.path.join(self._config.var_directory_path, "pid.lock")
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
        dburl = f"sqlite:///{self._config.var_directory_path}/foo.db"
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
        self._last_job_id_filename = os.path.join(self._config.var_directory_path, "last_job.txt")
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
                os.path.join(self._config.var_directory_path, "submission", f"{sid}.json")
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
                dbe=self._dbe,
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

    def _fetch_pending_job(self) -> Optional[ServerSideYSJSJob]:
        def _is_dependencies_cleared(_job: ServerSideYSJSJob) -> bool:
            with sessionmaker(bind=self._dbe)() as session:
                for depend_job in session.query(YSJSJobTable.job_id).filter(
                        YSJSJobTable.submission_id.in_(_job.submission.depends)
                ).all():
                    depend_job: YSJSJobTable
                    if depend_job.status != "finished":
                        return False
            return True

        if not self._job_queue_pending or self._state != "running":
            return None
        with self._job_queue_lock:
            if self._config.schedule_method == "FIFO":
                try:
                    smallest_jid = min(*self._job_queue_pending.keys())
                except TypeError:
                    return None
                job = self._job_queue_pending[smallest_jid]
                if (
                        _is_dependencies_cleared(job) and
                        job.submission.cpu < self._current_cpu and
                        job.submission.mem < self._current_mem
                ):
                    return self._job_queue_pending.pop(smallest_jid)
            elif self._config.schedule_method == "AGGRESSIVE":
                for i in list(self._job_queue_pending.keys()):
                    job = self._job_queue_pending[i]
                    if (
                            _is_dependencies_cleared(job) and
                            job.submission.cpu < self._current_cpu and
                            job.submission.mem < self._current_mem
                    ):
                        return self._job_queue_pending.pop(i)
        return None

    def job_cancel(self, job_id: int):
        """
        Cancel some job.


        :raise JobNotExistException: If the job does not exist in waiting queue
        """
        self._lh.debug("Cancel %d", job_id)
        if self._state == "starting":
            raise NotAvailableException
        with self._job_queue_lock:
            try:
                job_to_cancel = self._job_queue_pending.pop(job_id)
            except KeyError as e:
                raise JobNotExistException from e
        job_to_cancel.cancel()

    def job_send_signal(self, job_id: int, _signal: int):
        self._lh.debug("Kill -%d %d", _signal, job_id)
        if self._state == "starting":
            raise NotAvailableException
        try:
            self._job_queue_running[job_id].send_signal(_signal)
        except KeyError as e:
            raise JobNotExistException from e

    def job_kill(self, job_id: int):
        """
        Kill a running job.

        :raise JobNotExistException: If the job does not exist in running queue
        """
        self._lh.debug("Kill %d", job_id)
        if self._state == "starting":
            raise NotAvailableException
        try:
            self._job_queue_running[job_id].kill(self._config.kill_timeout)
        except KeyError as e:
            raise JobNotExistException from e

    def query(self, ) -> Iterable[int]:
        raise NotImplementedError

    def apply(
            self,
            job_ids: Iterable[int],
            operation: str,
            **extra_params
    ):
        job_ids = list(job_ids)
        if self._state == "starting":
            raise NotAvailableException
        if operation == "cancel":
            for job_id in job_ids:
                self.job_cancel(job_id)
        elif operation == "kill":
            for job_id in job_ids:
                self.job_kill(job_id)
        elif operation == "send_signal":
            for job_id in job_ids:
                self.job_send_signal(job_id, **extra_params)
        else:
            raise IllegalOperationException

    def run(self):
        self._lh.info("Started at http://localhost:%s", self._config.ysjs_port)
        self._state = "running"
        self._start_time = time.time()
        while self._state != "terminating":
            if len(self._job_queue_running) < self._config.max_concurrent_jobs:
                job_fetched = self._fetch_pending_job()
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

    @property
    def status(self) -> YSJSDStatus:
        return YSJSDStatus(
            state=self._state,
            current_cpu=self._current_cpu,
            current_mem=self._current_mem,
            pending_queue_length=len(self._job_queue_pending),
            running_queue_length=len(self._job_queue_running),
            uptime=time.time() - self._start_time
        )
