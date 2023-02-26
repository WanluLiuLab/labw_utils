from __future__ import annotations

from typing import Final, Union, Mapping, Any

AVAILABLE_SCHEDULING_METHOD = ("FIFO", "AGGRESSIVE")


class YSJSDConfig:
    _title: Final[str] = "ysjsd"
    _name: str
    _description: str
    _ysjs_port: str
    _var_directory_path: str
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
            var_directory_path: str,
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
        self._var_directory_path = var_directory_path
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
    def var_directory_path(self) -> str:
        return self._var_directory_path

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "name": self._name,
            "description": self._description,
            "ysjs_port": self._ysjs_port,
            "var_directory_path": self._var_directory_path,
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
