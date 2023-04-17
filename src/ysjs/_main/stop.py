import argparse

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List
from libysjs.operation import YSJSCluster

_lh = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--connection",
        required=False,
        help="YSJSD connection",
        nargs='?',
        type=str,
        action='store',
        default="http://localhost:8080"
    )
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    cl = YSJSCluster(conn=args.connection)
    _lh.info(
        "YSJS %s Cluster %s -- %s",
        cl.config.schedule_method, cl.config.name, cl.config.description
    )
    cl.stop()
    _lh.info(
        "Successfully stopped"
    )
