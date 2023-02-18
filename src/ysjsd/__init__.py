import json
import logging
import multiprocessing
import os.path
import platform
import sys
import threading
import time
import uuid
from typing import Mapping, Any, Final, List, Union, Optional, Dict

import flask
import psutil

import labw_utils

__version__ = "0.1.0"

from labw_utils.commonutils.io.safe_io import get_writer

from labw_utils.commonutils.stdlib_helper import logger_helper

from labw_utils.mlutils.io_helper import AbstractTOMLSerializable

AVAILABLE_SCHEDULING_METHOD=("FIFO", "AGGRESSIVE")

class YSJSDConfig(AbstractTOMLSerializable):
    _title: Final[str] = "ysjsd"
    _name: str
    _description: str
    _ysjs_port: str
    _var_directory: str
    _config_file_path: str
    _total_cpu: Union[int, float]
    _total_mem: Union[int, float]
    _schedule_method: str

    @staticmethod
    def _validate_versions(versions: Mapping[str, Any]) -> None:
        pass

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

    def __init__(
            self,
            name: str,
            description: str,
            ysjs_port: str,
            var_directory: str,
            config_file_path: str,
            total_cpu: Union[int, float],
            total_mem: Union[int, float],
            schedule_method: str
    ):
        self._name = name
        self._description = description
        self._ysjs_port = ysjs_port
        self._var_directory = var_directory
        self._config_file_path = config_file_path
        self._total_cpu = total_cpu
        self._total_mem = total_mem
        self._schedule_method = schedule_method

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
            "schedule_method": self._schedule_method
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
            schedule_method="FIFO"
        )

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

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


class YSJSSubmission:
    _job_id: str
    _job_name: str
    _job_description: str
    _queue_id: str
    _args: List[str]
    _cpu: Union[int, float]
    _mem: Union[int, float]
    _submission_time: float

    def __init__(
            self,
            job_id: str,
            job_name: str,
            job_description: str,
            queue_id: str,
            cpu: Union[int, float],
            mem: Union[int, float],
            args: List[str],
            submission_time: float
    ):
        self._job_id = job_id
        self._job_name = job_name
        self._cpu = cpu
        self._mem = mem
        self._args = args
        self._job_description = job_description
        self._queue_id = queue_id
        self._submission_time = submission_time

    @classmethod
    def new(
            cls,
            args: List[str],
            cpu: Union[int, float] = 1,
            mem: Union[int, float] = 1024 * 1024,
            job_name: str = "Unnamed",
            job_description: str = "No description",
            queue_id: str = "root"
    ):
        return cls(
            job_id=str(uuid.uuid4()),
            job_name=job_name,
            cpu=cpu,
            mem=mem,
            args=args,
            queue_id=queue_id,
            job_description=job_description,
            submission_time=time.time()
        )

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "job_id": self._job_id,
            "job_name": self._job_name,
            "cpu": self._cpu,
            "mem": self._mem,
            "args": self._args,
            "queue_id": self._queue_id,
            "job_description": self._job_description,
            "submission_time": self._submission_time
        }

    @property
    def submission_time(self) -> float:
        return self._submission_time

    @property
    def job_id(self) -> str:
        return self._job_id

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

    def to_json(self) -> str:
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, in_json: str):
        return cls.from_dict(json.loads(in_json))

    @property
    def cpu(self) -> Union[int, float]:
        return self._cpu

    @property
    def mem(self) -> Union[int, float]:
        return self._mem

    @property
    def queue_id(self) -> str:
        return self._queue_id


class YSJSDQueue:
    _queue_id: str
    _queue_name: str
    _queue_description: str
    _queue: List[YSJSSubmission]
    _lock: threading.Lock

    def __init__(
            self,
            queue_id: str,
            queue_name: str,
            queue_description: str
    ):
        self._lock = threading.Lock()
        self._queue_id = queue_id
        self._queue_name = queue_name
        self._queue_description = queue_description

    @classmethod
    def new(
            cls,
            queue_name: str = "Unnamed",
            queue_description: str = "No description"
    ):
        return cls(
            str(uuid.uuid4()),
            queue_name=queue_name,
            queue_description=queue_description
        )

    def add(self, submission: YSJSSubmission):
        with self._lock:
            self._queue.append(submission)

    def get_submission_aggressive(
            self,
            avail_cpu: Union[float, int],
            avail_mem: Union[float, int]
    ) -> Optional[YSJSSubmission]:
        with self._lock:
            for submission in self._queue:
                if submission.cpu < avail_cpu and submission.mem < avail_mem:
                    return submission
            return None

    def get_submission_fifo(
            self,
            avail_cpu: Union[float, int],
            avail_mem: Union[float, int]
    ) -> Optional[YSJSSubmission]:
        with self._lock:
            if not self._queue:
                return None
            else:
                submission = self._queue[0]
                if submission.cpu < avail_cpu and submission.mem < avail_mem:
                    self._queue.pop(0)
                    return submission

    @property
    def earliest_submission_time(self) -> Optional[float]:
        with self._lock:
            if not self._queue:
                return None
            else:
                return self._queue[0].submission_time

    @property
    def __len__(self):
        return len(self._queue)


class YSJSD(threading.Thread):
    _lock: threading.Lock
    _config: YSJSDConfig
    _current_cpu: Union[int, float]
    _current_mem: Union[int, float]
    _should_terminate: bool
    _schedule_method: str
    _lh: logging.Logger
    _queue: Dict[str, YSJSDQueue]

    def __init__(self, config: YSJSDConfig):
        super().__init__()
        self._lock = threading.Lock()
        self._config = config
        log_file_name = os.path.join(self._config.var_directory, "ysjsd.log")
        with get_writer(log_file_name) as w:
            pass  # Create file
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

    def add_job(self, submission: YSJSSubmission):
        if submission.queue_id in self._queue:
            self._queue[submission.queue_id].add(submission)
        else:
            pass

    def run(self):
        self._lh.info("Started")
        while not self._should_terminate:
            time.sleep(1)
        self._lh.info("Terminated")

    def terminate(self):
        self._lh.info("Received termination signal")
        self._should_terminate = True
