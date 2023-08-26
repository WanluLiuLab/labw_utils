import os
import tempfile

import pytest
from flask import Flask
from flask.testing import FlaskClient

from labw_utils.typing_importer import Tuple
from ysjsd.ds.ysjsd_config import ServerSideYSJSDConfig
from ysjsd.operation import YSJSD
from ysjsd.server import setup_globals


@pytest.fixture(scope="module", autouse=False)
def ysjsd_test_prep() -> Tuple[Flask, YSJSD, FlaskClient]:
    with tempfile.TemporaryDirectory() as test_dir:
        config_path = os.path.join(test_dir, "config.toml")
        var_path = os.path.join(test_dir, "var")
        config = ServerSideYSJSDConfig.new(config_file_path=config_path, var_directory_path=var_path)
        config.save(config_path)
        setup_globals(config)
        from ysjsd.server import global_flask_app
        from ysjsd.server import global_ysjsd

        yield global_flask_app, global_ysjsd, global_flask_app.test_client()
        global_ysjsd.terminate()
        global_ysjsd.join()
