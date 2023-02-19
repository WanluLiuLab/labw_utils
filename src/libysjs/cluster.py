from __future__ import annotations

import json
from typing import Iterable, Mapping, Union, Final, Any, List

import requests

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from libysjs.queue import YSJSDQueueConfig

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

    def job_submit(self):
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

    def queue_add(self):
        ...

    def queue_remove(self):
        ...

    def queue_pause(self):
        ...

    def queue_resume(self):
        ...

    def queue_kill(self):
        ...

    def queue_cancel(self):
        ...

    @property
    def queue_names(self) -> Iterable[str]:
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
    _queue_conf: List[YSJSDQueueConfig]

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
            queue_conf: List[YSJSDQueueConfig]
    ):
        self._name = name
        self._description = description
        self._ysjs_port = ysjs_port
        self._var_directory = var_directory
        self._config_file_path = config_file_path
        self._total_cpu = total_cpu
        self._total_mem = total_mem
        self._schedule_method = schedule_method
        self._queue_conf = queue_conf

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
            "queue_conf": list(map(YSJSDQueueConfig.to_dict, self._queue_conf))
        }

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        in_dict = dict(in_dict)
        queue_conf = list(map(lambda _in_dict: YSJSDQueueConfig.from_dict(**_in_dict), in_dict.pop("queue_conf")))
        return cls(**in_dict, queue_conf=queue_conf)

