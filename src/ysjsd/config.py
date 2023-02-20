from __future__ import annotations

import multiprocessing
import os
import platform
import sys
import time
from typing import Mapping, Any

import flask
import psutil
import sqlalchemy

import labw_utils
from labw_utils.commonutils.serializer.toml import AbstractTOMLSerializable
from labw_utils.commonutils.stdlib_helper import logger_helper
from libysjs import __version__
from libysjs.cluster import CoreYSJSDConfig, AVAILABLE_SCHEDULING_METHOD


class ServerSideYSJSDConfig(CoreYSJSDConfig, AbstractTOMLSerializable):
    @staticmethod
    def _validate_versions(versions: Mapping[str, Any]) -> None:
        _lh = logger_helper.get_logger("YSJSD - LOADING CONFIG")
        compile_time_version_dict = dict(versions)
        run_time_version_dict = ServerSideYSJSDConfig._dump_versions()
        for version_key, run_time_version_value in run_time_version_dict.items():
            try:
                compile_time_version_value = compile_time_version_dict.pop(version_key)
            except KeyError:
                compile_time_version_value = None
            if compile_time_version_value != run_time_version_value:
                _lh.warning(
                    "Package %s have different version information: Compile (%s) != Run (%s)",
                    version_key, compile_time_version_value, run_time_version_value
                )
        for remaining_compile_time_version_key, remaining_compile_time_version_value in compile_time_version_dict.items():
            _lh.warning(
                "Package %s have different version information: Compile (%s) != Run (%s)",
                remaining_compile_time_version_key, remaining_compile_time_version_value, None
            )

    @staticmethod
    def _dump_versions() -> Mapping[str, Any]:
        return {
            "labw_utils": labw_utils.__version__,
            "ysjsd": __version__,
            "flask": flask.__version__,
            "sqlalchemy": sqlalchemy.__version__,
            "python": ".".join(map(str, (sys.version_info.major, sys.version_info.minor, sys.version_info.micro)))
        }

    @staticmethod
    def dump_versions() -> Mapping[str, Any]:
        return ServerSideYSJSDConfig._dump_versions()

    @staticmethod
    def _dump_metadata() -> Mapping[str, Any]:
        return {
            "time_gmt": time.asctime(time.gmtime()),
            "time_local": time.asctime(time.localtime()),
            "platform": platform.uname()._asdict(),
            "interpreter": {
                "executable": sys.executable,
                "path": sys.path,
                "implementation": sys.implementation.name
            },
            "env": dict(os.environ),
            "cwd": os.getcwd()
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
            schedule_method="FIFO",
            max_concurrent_jobs=4096,
            kill_timeout=3
        )

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
