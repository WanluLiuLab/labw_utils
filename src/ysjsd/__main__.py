import os
import argparse
import sys

from typing import Optional, List

import flask

from ysjsd import YSJSDConfig

app_dir = os.path.dirname(os.path.abspath(__file__))
global_config: Optional[YSJSDConfig] = None
app = flask.Flask(
    "ylsjsd",
    template_folder=os.path.join(app_dir, "templates")
)


def get_or_create_db(db_filepath: str):
    ...


def create_rest_service(ylsjs_port: str):
    ...


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--generate_default_config',
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


@app.route('/ylsjsd/api/v1.0/config', methods=['GET'])
def serve_config():
    if global_config is None:
        raise ValueError
    return flask.jsonify(**global_config.to_dict())


@app.route('/', methods=['GET'])
def serve_frontend():
    return flask.render_template("frontpage.html")


def start(config: YSJSDConfig):
    global global_config
    global_config = config
    app.run(
        host="0.0.0.0",
        port=int(config.ysjs_port)
    )


if __name__ == "__main__":
    args = _parse_args(sys.argv[1:])
    if args.generate_default_config:
        YSJSDConfig(
            "ylsjs_cluster",
            "No description",
            "8080",
            os.path.join(os.path.abspath("."), "var"),
            os.path.abspath(args.config)
        ).to_config(args.config)
        sys.exit(0)
    start(YSJSDConfig.from_config(args.config))
