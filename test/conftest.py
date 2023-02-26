"""
test_tetgs.py -- Unit test of corresponding module.
"""
import os
import tempfile
from typing import Tuple

import pytest
from flask import Flask

from labw_utils.commonutils import shell_utils
from labw_utils.commonutils.stdlib_helper import logger_helper
from libysjs.ds.ysjsd_config import YSJSDConfig
from ysjsd import APP_NAME, APP_DIR
from ysjsd.ds.ysjsd_config import ServerSideYSJSDConfig
from ysjsd.operation import YSJSD


class SessionTestInfo:
    __slots__ = (
        "base_test_dir"
    )
    base_test_dir: str

    def __init__(self):
        self.base_test_dir = tempfile.mkdtemp()
        lh = logger_helper.get_logger("TEST_TETGS")
        lh.info(f"Test dir {self.base_test_dir}")

    def teardown(self):
        shell_utils.rm_rf(self.base_test_dir)


class ModuleTestInfo:
    __slots__ = (
        "name",
        "path"
    )
    path: str

    def __init__(self, base_test_dir: str, name: str):
        self.name = name
        self.path = os.path.join(base_test_dir, name)
        shell_utils.mkdir_p(self.path)

    def teardown(self):
        shell_utils.rm_rf(self.path)


@pytest.fixture(scope="session")
def initialize_session():
    session_test_info = SessionTestInfo()
    yield session_test_info
    session_test_info.teardown()


@pytest.fixture(scope="module", autouse=False)
def ysjsd_test_prep() -> Tuple[Flask, YSJSD, YSJSDConfig]:
    with tempfile.TemporaryDirectory() as test_dir:
        config_path = os.path.join(test_dir, "config.toml")
        var_path = os.path.join(test_dir, "var")
        config = ServerSideYSJSDConfig.new(config_file_path=config_path, var_directory_path=var_path)
        config.save(config_path)
        app = Flask(
            APP_NAME,
            template_folder=os.path.join(APP_DIR, "templates")
        )
        ysjsd = YSJSD(config)
        yield app, ysjsd, config
        # Destroy

@pytest.fixture()
def ysjsd_client(ysjsd_test_prep):
    app, _, _ = ysjsd_test_prep
    return app.test_client()
