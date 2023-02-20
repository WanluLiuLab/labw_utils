from __future__ import annotations

import json
from typing import Mapping, Union, Final, Any

import requests

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from libysjs.submission import YSJSSubmission

_lh = get_logger(__name__)

AVAILABLE_SCHEDULING_METHOD = ("FIFO", "AGGRESSIVE")


class YSJSCluster:
    _config: CoreYSJSDConfig
    _conn: str

    def __init__(self, conn: str):
        resp = requests.get(f"{conn}/ysjsd/api/v1.0/config")
        if resp.status_code != 200:
            exit(1)
        else:
            try:
                self._config = CoreYSJSDConfig.from_dict(json.loads(resp.text))
            except Exception:
                exit(1)
        _lh.debug("Successfully GET YSJSD configuration")
        self._conn = conn

    def submit(self, submission: YSJSSubmission):
        resp = requests.post(
            f"{self._conn}/ysjsd/api/v1.0/submit",
            data=json.dumps(submission.to_dict())
        )
        if resp.status_code != 200:
            exit(1)
        _lh.debug("Successfully POST submission %s", submission.submission_name)

    def stop(self):
        ...

    def job_send_signal(self):
        ...

    def job_kill(self):
        ...

    def job_pause(self):
        ...

    def job_resume(self):
        ...

    def job_cancel(self):
        ...

    @property
    def cluster_load(self) -> YSJSDLoad:
        resp = requests.get(f"{self._conn}/ysjsd/api/v1.0/load")
        if resp.status_code != 200:
            exit(1)
        else:
            try:
                load = YSJSDLoad.from_dict(json.loads(resp.text))
            except Exception:
                exit(1)
        _lh.debug("Successfully GET YSJSD load")
        return load

    @property
    def config(self) -> CoreYSJSDConfig:
        return self._config


class YSJSDLoad:
    _real_avail_cpu: float
    _real_total_cpu: int
    _real_avail_mem: int
    _real_total_mem: int

    def __init__(
            self,
            real_avail_cpu: float,
            real_total_cpu: int,
            real_avail_mem: int,
            real_total_mem: int
    ):
        self._real_avail_cpu = real_avail_cpu
        self._real_total_mem = real_total_mem
        self._real_avail_mem = real_avail_mem
        self._real_total_cpu = real_total_cpu

    def to_dict(self) -> Mapping[str, Union[int, float]]:
        return {
            "real_avail_cpu": self._real_avail_cpu,
            "real_total_cpu": self._real_total_cpu,
            "real_avail_mem": self._real_avail_mem,
            "real_total_mem": self._real_total_mem
        }

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Union[int, float]]):
        return cls(**in_dict)

    @property
    def real_total_cpu(self) -> int:
        return self._real_total_cpu

    @property
    def real_total_mem(self) -> int:
        return self._real_total_mem

    @property
    def real_avail_mem(self) -> int:
        return self._real_avail_mem

    @property
    def real_avail_cpu(self) -> float:
        return self._real_avail_cpu


class CoreYSJSDConfig:
    _title: Final[str] = "ysjsd"
    _name: str
    _description: str
    _ysjs_port: str
    _var_directory: str
    _config_file_path: str
    _total_cpu: Union[int, float]
    _total_mem: Union[int, float]
    _schedule_method: str
    _max_concurrent_jobs: int
    _kill_timeout: float

    def __init__(
            self,
            name: str,
            description: str,
            ysjs_port: str,
            var_directory: str,
            config_file_path: str,
            total_cpu: Union[int, float],
            total_mem: Union[int, float],
            schedule_method: str,
            max_concurrent_jobs: int,
            kill_timeout: float
        ):
        self._name = name
        self._description = description
        self._ysjs_port = ysjs_port
        self._var_directory = var_directory
        self._config_file_path = config_file_path
        self._total_cpu = total_cpu
        self._total_mem = total_mem
        self._schedule_method = schedule_method
        self._max_concurrent_jobs = max_concurrent_jobs
        self._kill_timeout = kill_timeout


    @property
    def max_concurrent_jobs(self) -> int:
        return self._max_concurrent_jobs
    @property
    def kill_timeout(self) -> float:
        return self._kill_timeout
    @property
    def name(self) -> str:
        return self._name

    @property
    def description(self) -> str:
        return self._description

    @property
    def ysjs_port(self) -> str:
        return self._ysjs_port

    @property
    def total_cpu(self) -> Union[int, float]:
        return self._total_cpu

    @property
    def total_mem(self) -> Union[int, float]:
        return self._total_mem

    @property
    def schedule_method(self) -> str:
        return self._schedule_method

    @property
    def var_directory(self) -> str:
        return self._var_directory

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "name": self._name,
            "description": self._description,
            "ysjs_port": self._ysjs_port,
            "var_directory": self._var_directory,
            "config_file_path": self._config_file_path,
            "total_cpu": self._total_cpu,
            "total_mem": self._total_mem,
            "schedule_method": self._schedule_method,
            "max_concurrent_jobs": self._max_concurrent_jobs,
            "kill_timeout": self._kill_timeout
        }

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

