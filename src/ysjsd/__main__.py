import argparse
import os
import sys
from typing import List

from labw_utils.commonutils.stdlib_helper import logger_helper
from ysjsd.ds.ysjsd_config import ServerSideYSJSDConfig
from ysjsd.server import start

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


if __name__ == "__main__":
    args = _parse_args(sys.argv[1:])
    if args.generate_default_config:
        ServerSideYSJSDConfig.new(args.config).save(args.config)
        _lh.info("Configure file generated at %s", args.config)
        sys.exit(0)
    start(ServerSideYSJSDConfig.load(args.config))
