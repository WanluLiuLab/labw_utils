from __future__ import annotations

import json

from labw_utils import UnmetDependenciesError
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Mapping, Union
from libysjs.ds.ysjs_submission import YSJSSubmission
from libysjs.ds.ysjsd_config import YSJSDConfig
from libysjs.ds.ysjsd_status import YSJSDStatus

try:
    import pytest

    requests = pytest.importorskip("requests")
except ImportError:
    pytest = None
    try:
        import requests
    except ImportError:
        raise UnmetDependenciesError("requests")

_lh = get_logger(__name__)


class YSJSException(RuntimeError):
    ...


class ClusterNotUpException(YSJSException):
    ...


class MalformedResponseException(YSJSException):
    ...


class YSJSCluster:
    _config: YSJSDConfig
    _conn: str

    def __init__(self, conn: str):
        try:
            resp = requests.get(f"{conn}/ysjsd/api/v1.0/config")
        except requests.RequestException as e:
            raise ClusterNotUpException from e
        if resp.status_code != 200:
            raise MalformedResponseException
        else:
            try:
                self._config = YSJSDConfig.from_dict(json.loads(resp.text))
            except json.JSONDecodeError as e:
                raise MalformedResponseException from e
        _lh.debug("Successfully GET YSJSD configuration")
        self._conn = conn

    def submit(self, submission: YSJSSubmission) -> int:
        try:
            resp = requests.post(f"{self._conn}/ysjsd/api/v1.0/submit", data=json.dumps(submission.to_dict()))
        except requests.RequestException as e:
            raise ClusterNotUpException from e
        if resp.status_code != 200:
            raise MalformedResponseException
        try:
            job_id = int(resp.text)
        except ValueError as e:
            raise MalformedResponseException from e
        _lh.debug("Successfully POST submission %s -> job (%d)", submission.submission_name, job_id)
        return job_id

    def stop(self):
        try:
            resp = requests.post(f"{self._conn}/ysjsd/api/v1.0/stop")
        except requests.RequestException as e:
            raise ClusterNotUpException from e
        if resp.status_code != 200:
            raise MalformedResponseException
        _lh.debug("Successfully POST stop")

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
        try:
            resp = requests.get(f"{self._conn}/ysjsd/api/v1.0/load")
        except requests.RequestException as e:
            raise ClusterNotUpException from e
        if resp.status_code != 200:
            raise MalformedResponseException
        else:
            try:
                load = YSJSDLoad.from_dict(json.loads(resp.text))
            except json.JSONDecodeError as e:
                raise MalformedResponseException from e
        _lh.debug("Successfully GET YSJSD load")
        return load

    @property
    def cluster_status(self) -> YSJSDStatus:
        try:
            resp = requests.get(f"{self._conn}/ysjsd/api/v1.0/status")
        except requests.RequestException as e:
            raise ClusterNotUpException from e
        if resp.status_code != 200:
            raise MalformedResponseException
        else:
            try:
                load = YSJSDStatus.from_dict(json.loads(resp.text))
            except json.JSONDecodeError as e:
                raise MalformedResponseException from e
        _lh.debug("Successfully GET YSJSD status")
        return load

    @property
    def config(self) -> YSJSDConfig:
        return self._config


class YSJSDLoad:
    _real_avail_cpu: float
    _real_total_cpu: int
    _real_avail_mem: int
    _real_total_mem: int

    def __init__(self, real_avail_cpu: float, real_total_cpu: int, real_avail_mem: int, real_total_mem: int):
        self._real_avail_cpu = real_avail_cpu
        self._real_total_mem = real_total_mem
        self._real_avail_mem = real_avail_mem
        self._real_total_cpu = real_total_cpu

    def to_dict(self) -> Mapping[str, Union[int, float]]:
        return {
            "real_avail_cpu": self._real_avail_cpu,
            "real_total_cpu": self._real_total_cpu,
            "real_avail_mem": self._real_avail_mem,
            "real_total_mem": self._real_total_mem,
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
