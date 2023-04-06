import base64
import json
import logging
import os
import signal
from typing import Optional, Tuple, Union

import flask
from gevent import pywsgi

from labw_utils.commonutils.stdlib_helper import logger_helper
from libysjs.ds.ysjs_submission import YSJSSubmission
from ysjsd.ds.ysjsd_config import ServerSideYSJSDConfig
from ysjsd.operation import YSJSD, JobNotExistException

# Global Constants
APP_DIR = os.path.dirname(os.path.abspath(__file__))
APP_NAME = "YSJSD BACKEND"
NOT_READY = ("Not ready\n", 500)

# Global Variables
global_config: Optional[ServerSideYSJSDConfig] = None
global_ysjsd: Optional[YSJSD] = None
global_server: Optional[pywsgi.WSGIServer] = None
global_flask_app = flask.Flask(
    APP_NAME,
    template_folder=os.path.join(APP_DIR, "templates")
)

# Typings
ResponseType = Tuple[Union[str, flask.Response], int]

# Create Logger
stream_handler = logging.StreamHandler()
stream_handler.setLevel('INFO')
stream_handler.setFormatter(logger_helper.get_formatter(stream_handler.level))

logging.basicConfig(
    handlers=[
        stream_handler
    ],
    force=True,
    level=logger_helper.TRACE
)

_lh = logger_helper.get_logger("YSJSD BACKEND")


@global_flask_app.route('/ysjsd/api/v1.0/config', methods=['GET'])
def serve_config() -> ResponseType:
    global global_config
    if global_config is None:
        return NOT_READY
    return flask.jsonify(**global_config.to_dict()), 200


@global_flask_app.route('/ysjsd/api/v1.0/load', methods=['GET'])
def serve_load() -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    return flask.jsonify(**global_ysjsd.real_load.to_dict()), 200


@global_flask_app.route('/ysjsd/api/v1.0/status', methods=['GET'])
def serve_status() -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    return flask.jsonify(**global_ysjsd.status.to_dict()), 200


@global_flask_app.route('/ysjsd/api/v1.0/submission/<submission_id>', methods=['GET'])
def serve_submission(submission_id: str) -> ResponseType:
    ...


@global_flask_app.route('/ysjsd/api/v1.0/job/<int:job_id>', methods=['GET'])
def serve_job(job_id: int) -> ResponseType:
    ...


@global_flask_app.route('/ysjsd/api/v1.0/job/<int:job_id>/cancel', methods=['POST'])
def cancel(job_id: int) -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    try:
        global_ysjsd.job_cancel(job_id)
        return f"Cancel {job_id}\n", 200
    except JobNotExistException:
        return f"Cancel {job_id} Failure -- Job not exist\n", 500


@global_flask_app.route('/ysjsd/api/v1.0/job/<int:job_id>/send_signal/<int:_signal>', methods=['POST'])
def send_signal(job_id: int, _signal: int) -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    try:
        global_ysjsd.job_send_signal(job_id, _signal)
        return f"Send signal {_signal} to {job_id}\n", 200
    except JobNotExistException:
        return f"Send signal {_signal} to {job_id} Failure -- Job not exist\n", 500


@global_flask_app.route('/ysjsd/api/v1.0/job/<int:job_id>/kill', methods=['POST'])
def kill(job_id: int) -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    try:
        global_ysjsd.job_kill(job_id)
        return f"Kill {job_id}\n", 200
    except JobNotExistException:
        return f"Kill {job_id} Failure -- Job not exist\n", 500


@global_flask_app.route('/ysjsd/api/v1.0/stop', methods=['POST'])
def stop() -> ResponseType:
    global global_server, global_ysjsd
    global_server.stop()
    global_ysjsd.terminate()
    global_ysjsd.join()
    return "STOPPED/n", 200


@global_flask_app.route('/ysjsd/api/v1.0/submit', methods=['POST'])
def receive_submission() -> ResponseType:
    global global_ysjsd
    data = flask.request.get_data()
    try:
        submission = YSJSSubmission.from_dict(
            json.loads(str(data, encoding="UTF8"))
        )
    except Exception as e:
        err_message = f"{str(e)} when parse submission {str(base64.b64encode(data), encoding='UTF8')}"
        return err_message, 500
    try:
        ret_jid = global_ysjsd.receive_submission(submission)
        return str(ret_jid), 200
    except ValueError as e:
        err_message = f"{str(e)} when parse submission {str(base64.b64encode(data), encoding='UTF8')}"
        return err_message, 500


@global_flask_app.route('/', methods=['GET'])
def serve_frontend() -> ResponseType:
    return flask.render_template("frontpage.html"), 200


def setup_globals(config: ServerSideYSJSDConfig):
    global global_config, global_ysjsd, global_server
    global_config = config
    global_ysjsd = YSJSD(global_config)
    global_ysjsd.start()
    global_flask_app.logger.handlers.clear()
    global_flask_app.logger.setLevel(logger_helper.TRACE)
    frontend_logger_file_handler = logging.FileHandler(
        os.path.join(global_config.var_directory_path, "ysjsd_pywsgi.log")
    )
    frontend_logger_file_handler.setLevel(logger_helper.TRACE)
    frontend_logger_file_handler.setFormatter(logger_helper.get_formatter(frontend_logger_file_handler.level))
    global_flask_app.logger.addHandler(frontend_logger_file_handler)
    signal.signal(signal.SIGINT, lambda x, y: stop())
    signal.signal(signal.SIGTERM, lambda x, y: stop())
    signal.signal(signal.SIGHUP, lambda x, y: stop())
    global_server = pywsgi.WSGIServer(
        ("0.0.0.0", int(global_config.ysjs_port)),
        application=global_flask_app,
        log=pywsgi.LoggingLogAdapter(global_flask_app.logger, level=logging.DEBUG),
        error_log=None
    )


def start(config: ServerSideYSJSDConfig):
    global global_config, global_ysjsd, global_server
    setup_globals(config)
    global_server.serve_forever()
