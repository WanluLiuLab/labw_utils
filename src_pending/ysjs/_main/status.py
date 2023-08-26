import argparse

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List
from libysjs.operation import YSJSCluster
from libysjs.utils import scale_si

_lh = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--connection",
        required=False,
        help="YSJSD connection",
        nargs="?",
        type=str,
        action="store",
        default="http://localhost:8080",
    )
    parser.add_argument("--show_conf", required=False, help="Show cluster configuration", action="store_true")
    parser.add_argument("--show_load", required=False, help="Show real cluster load", action="store_true")
    parser.add_argument("--show_status", required=False, help="Show real cluster status", action="store_true")
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    cl = YSJSCluster(conn=args.connection)
    _lh.info("YSJS %s Cluster %s -- %s", cl.config.schedule_method, cl.config.name, cl.config.description)
    if args.show_conf:
        total_mem_si, total_mem_si_prefix = scale_si(cl.config.total_mem)
        _lh.info(
            "Configured Total Resources: CPU %.2f, Memory %.2f %sB (%.2f)",
            cl.config.total_cpu,
            total_mem_si,
            total_mem_si_prefix,
            cl.config.total_mem,
        )
    if args.show_load:
        current_load = cl.cluster_load
        real_total_mem_si, real_total_mem_si_prefix = scale_si(current_load.real_total_mem)
        _lh.info(
            "Real Total Resources: CPU %.2f, Memory %.2f %sB (%.2f)",
            current_load.real_total_cpu,
            real_total_mem_si,
            real_total_mem_si_prefix,
            current_load.real_total_mem,
        )
        real_avail_mem_si, real_avail_mem_prefix = scale_si(current_load.real_avail_mem)
        _lh.info(
            "Real Available Resources: CPU %.2f, Memory %.2f %sB (%.2f)",
            current_load.real_avail_cpu,
            real_avail_mem_si,
            real_avail_mem_prefix,
            current_load.real_avail_mem,
        )
    if args.show_status:
        current_status = cl.cluster_status
