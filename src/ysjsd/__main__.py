import argparse
import base64
import json
import logging
import os
import signal
import sys
from typing import Optional, List, Tuple, Union

import flask
from gevent import pywsgi
from gevent.pywsgi import WSGIServer

from labw_utils.commonutils.stdlib_helper import logger_helper
from libysjs.ds.ysjs_submission import YSJSSubmission
from ysjsd import APP_DIR, APP_NAME, NOT_READY
from ysjsd.ds.ysjsd_config import ServerSideYSJSDConfig
from ysjsd.operation import YSJSD, JobNotExistException

global_config: Optional[ServerSideYSJSDConfig] = None
global_ysjsd: Optional[YSJSD] = None
global_server: WSGIServer

app = flask.Flask(
    APP_NAME,
    template_folder=os.path.join(APP_DIR, "templates")
)

stream_handler = logging.StreamHandler()
stream_handler.setLevel('INFO')
stream_handler.setFormatter(logger_helper.get_formatter(stream_handler.level))

ResponseType = Tuple[Union[str, flask.Response], int]

logging.basicConfig(
    handlers=[
        stream_handler
    ],
    force=True,
    level=logger_helper.TRACE
)
_lh = logger_helper.get_logger("YSJSD BACKEND")


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g", '--generate_default_config',
        required=False, help="Generate default config and exit",
        action="store_true"
    )
    parser.add_argument(
        '-e', '--config', required=False,
        help="Config file path", nargs='?',
        type=str, action='store',
        default=os.path.join(os.path.abspath("."), "config.toml")
    )
    return parser.parse_args(args)


@app.route('/ysjsd/api/v1.0/config', methods=['GET'])
def serve_config() -> ResponseType:
    global global_config
    if global_config is None:
        return NOT_READY
    return flask.jsonify(**global_config.to_dict()), 200


@app.route('/ysjsd/api/v1.0/load', methods=['GET'])
def serve_load() -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    return flask.jsonify(**global_ysjsd.real_load.to_dict()), 200


@app.route('/ysjsd/api/v1.0/status', methods=['GET'])
def serve_status() -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    return flask.jsonify(**global_ysjsd.status.to_dict()), 200


@app.route('/ysjsd/api/v1.0/submission/<submission_id>', methods=['GET'])
def serve_submission(submission_id: str) -> ResponseType:
    ...


@app.route('/ysjsd/api/v1.0/job/<int:job_id>', methods=['GET'])
def serve_job(job_id: int) -> ResponseType:
    ...


@app.route('/ysjsd/api/v1.0/job/<int:job_id>/cancel', methods=['POST'])
def cancel(job_id: int) -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    try:
        global_ysjsd.job_cancel(job_id)
        return f"Cancel {job_id}\n", 200
    except JobNotExistException:
        return f"Cancel {job_id} Failure -- Job not exist\n", 500


@app.route('/ysjsd/api/v1.0/job/<int:job_id>/send_signal/<int:_signal>', methods=['POST'])
def send_signal(job_id: int, _signal: int) -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    try:
        global_ysjsd.job_send_signal(job_id, _signal)
        return f"Send signal {_signal} to {job_id}\n", 200
    except JobNotExistException:
        return f"Send signal {_signal} to {job_id} Failure -- Job not exist\n", 500


@app.route('/ysjsd/api/v1.0/job/<int:job_id>/kill', methods=['POST'])
def kill(job_id: int) -> ResponseType:
    global global_ysjsd
    if global_ysjsd is None:
        return NOT_READY
    try:
        global_ysjsd.job_kill(job_id)
        return f"Kill {job_id}\n", 200
    except JobNotExistException:
        return f"Kill {job_id} Failure -- Job not exist\n", 500


@app.route('/ysjsd/api/v1.0/stop', methods=['POST'])
def stop() -> ResponseType:
    global global_server, global_ysjsd
    global_server.stop()
    global_ysjsd.terminate()
    global_ysjsd.join()
    return "STOPPED/n", 200


@app.route('/ysjsd/api/v1.0/submit', methods=['POST'])
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


@app.route('/', methods=['GET'])
def serve_frontend() -> ResponseType:
    return flask.render_template("frontpage.html"), 200


def start(config: ServerSideYSJSDConfig):
    global global_config, global_ysjsd, global_server
    global_config = config
    global_ysjsd = YSJSD(global_config)
    global_ysjsd.start()
    # Update the logger
    app.logger.handlers.clear()
    app.logger.setLevel(logger_helper.TRACE)
    frontend_logger_file_handler = logging.FileHandler(
        os.path.join(global_config.var_directory_path, "ysjsd_pywsgi.log")
    )
    frontend_logger_file_handler.setLevel(logger_helper.TRACE)
    frontend_logger_file_handler.setFormatter(logger_helper.get_formatter(frontend_logger_file_handler.level))
    app.logger.addHandler(frontend_logger_file_handler)
    signal.signal(signal.SIGINT, lambda x, y: stop())
    signal.signal(signal.SIGTERM, lambda x, y: stop())
    signal.signal(signal.SIGHUP, lambda x, y: stop())
    global_server = pywsgi.WSGIServer(
        ("0.0.0.0", int(global_config.ysjs_port)),
        application=app,
        log=pywsgi.LoggingLogAdapter(app.logger, level=logging.DEBUG),
        error_log=None
    )
    global_server.serve_forever()


if __name__ == "__main__":
    args = _parse_args(sys.argv[1:])
    if args.generate_default_config:
        ServerSideYSJSDConfig.new(args.config).save(args.config)
        _lh.info("Configure file generated at %s", args.config)
        sys.exit(0)
    start(ServerSideYSJSDConfig.load(args.config))
