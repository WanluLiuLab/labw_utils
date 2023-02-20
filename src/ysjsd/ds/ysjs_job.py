from __future__ import annotations

import signal
import subprocess
import threading
import time
from typing import Optional

import psutil

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from libysjs.ds.ysjs_job import YSJSJob

_lh = get_logger("YSJSD")


class ServerSideYSJSJob(threading.Thread, YSJSJob):
    _p: Optional[subprocess.Popen]

    def __init__(self, **kwargs):
        threading.Thread.__init__(self)
        YSJSJob.__init__(self, **kwargs)
        self._p = None

    def run(self):
        _lh.info("Job %d starting", self._job_id)
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
        self._pid = self._p.pid
        self._start_time = time.time()
        _lh.info("Job %d running", self._job_id)
        self._status = "running"
        self._retv = self._p.wait()
        self._status = "finished"
        _lh.info("Job %d finished with exit value %d", self._job_id, self._retv)
        self._terminate_time = time.time()

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

    def kill(self, timeout: float):
        """
        Recursively terminate a process tree
        """

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
