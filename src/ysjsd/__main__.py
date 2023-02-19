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
from libysjs.submission import YSJSSubmission
from ysjsd import YSJSD, ServerSideYSJSDConfig

APP_DIR = os.path.dirname(os.path.abspath(__file__))
APP_NAME = "YSJSD BACKEND"

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
    if global_config is None:
        raise ValueError
    return flask.jsonify(**global_config.to_dict()), 200


@app.route('/ysjsd/api/v1.0/load', methods=['GET'])
def serve_load() -> ResponseType:
    if global_ysjsd is None:
        raise ValueError
    return flask.jsonify(**global_ysjsd.get_real_load().to_dict()), 200


@app.route('/ysjsd/api/v1.0/queue_names', methods=['GET'])
def serve_queue_names() -> ResponseType:
    ...


@app.route('/ysjsd/api/v1.0/queue/<queue_name>', methods=['GET'])
def serve_queue(queue_name: str) -> ResponseType:
    ...


@app.route('/ysjsd/api/v1.0/job_ids', methods=['GET'])
def serve_job_ids() -> ResponseType:
    ...


@app.route('/ysjsd/api/v1.0/submission/<job_id>', methods=['GET'])
def serve_submission(job_id: str) -> ResponseType:
    ...



@app.route('/ysjsd/api/v1.0/stop', methods=['POST'])
def stop() -> ResponseType:
    global global_server, global_ysjsd
    global_server.stop()
    global_ysjsd.terminate()
    global_ysjsd.join()
    return "STOPPED/n", 200


@app.route('/ysjsd/api/v1.0/job/submit', methods=['POST'])
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
        global_ysjsd.receive_submission(submission)
        return f"added job_id {submission.job_id}", 200
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
    app.logger.setLevel(logger_helper.TRACE)
    frontend_logger_file_handler = logging.FileHandler(
        os.path.join(global_config.var_directory, "ysjsd_backend.log")
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
        log=pywsgi.LoggingLogAdapter(app.logger)
    )
    global_server.serve_forever()


if __name__ == "__main__":
    args = _parse_args(sys.argv[1:])
    if args.generate_default_config:
        ServerSideYSJSDConfig.new(
            args.config
        ).save(args.config)
        _lh.info("Configure file generated at %s", args.config)
        sys.exit(0)
    start(ServerSideYSJSDConfig.load(args.config))
