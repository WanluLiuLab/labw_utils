from __future__ import annotations

import signal
import subprocess
import threading
import time

import psutil
import sqlalchemy.engine
from sqlalchemy.orm import sessionmaker

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Optional, Any, Mapping
from libysjs.ds.ysjs_job import YSJSJob
from ysjsd.orm.ysjs_job_table import YSJSJobTable

_lh = get_logger("YSJSD")


class ServerSideYSJSJob(threading.Thread, YSJSJob):
    _p: Optional[subprocess.Popen]
    _dbe: sqlalchemy.engine.Engine
    _db_write_lock: threading.Lock

    def __init__(self, dbe: sqlalchemy.engine.Engine, db_write_lock: threading.Lock, **kwargs):
        threading.Thread.__init__(self)
        YSJSJob.__init__(self, **kwargs)
        self._p = None
        self._dbe = dbe
        self._db_write_lock = db_write_lock
        with self._db_write_lock:
            with sessionmaker(bind=self._dbe)() as session:
                session.add(YSJSJobTable.from_job(self))
                session.commit()

    def _db_update(self, update_dict: Mapping[str, Any]):
        with sessionmaker(bind=self._dbe)() as session:
            session.query(YSJSJobTable).filter(YSJSJobTable.job_id == self.job_id).update(update_dict)
            session.commit()

    @classmethod
    def new(cls, dbe: sqlalchemy.engine.Engine, db_write_lock: threading.Lock, **kwargs):
        return cls(dbe=dbe, db_write_lock=db_write_lock, **YSJSJob.new(**kwargs).to_dict())

    def run(self):
        _lh.info("Job %d starting", self._job_id)
        self._db_update({"status": self._status})
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
            args=[self._submission.shell_path, self._submission.script_path],
            cwd=self._submission.cwd,
            env=self._submission.env,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            close_fds=True,
        )
        self._pid = self._p.pid
        self._start_time = time.time()
        _lh.info("Job %d running", self._job_id)
        self._status = "running"
        self._db_update({"status": self._status, "pid": self._pid, "start_time": self.start_time})
        self._retv = self._p.wait()
        self._status = "finished"
        self._terminate_time = time.time()
        self._db_update({"status": self._status, "retv": self._retv, "terminate_time": self._terminate_time})
        _lh.info("Job %d finished with exit value %d", self._job_id, self._retv)

    def send_signal(self, _signal: int):
        if self._p is None:
            raise  # TODO
        self._p.send_signal(_signal)

    @property
    def pid(self) -> Optional[int]:
        if self._p is not None:
            return self._p.pid
        else:
            return None

    def cancel(self):
        self._status = "canceled"
        self._db_update({"status": self._status})

    def kill(self, timeout: float):
        """
        Recursively terminate a process tree
        """
        if self._p is None:
            raise  # TODO

        def _kill(_signal: int):
            try:
                p = psutil.Process(self._p.pid)
                children = p.children(recursive=True)
            except psutil.Error:
                return
            processes_to_be_killed = [p, *children]
            for process_to_be_killed in processes_to_be_killed:
                try:
                    process_to_be_killed.send_signal(_signal)
                except psutil.Error:
                    pass

        _kill(signal.SIGTERM)
        time.sleep(timeout)
        _kill(signal.SIGKILL)
