import argparse
import datetime
import os

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List, Union
from libysjs.ds.ysjs_submission import (
    YSJSSubmission,
    DEFAULT_SUBMISSION_NAME,
    DEFAULT_SUBMISSION_DESCRIPTION,
    DEFAULT_SUBMISSION_CPU,
    DEFAULT_SUBMISSION_MEM,
)
from libysjs.operation import YSJSCluster

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
    parser.add_argument(
        "--script_path", required=True, help="Path of script to be executed", nargs="?", type=str, action="store"
    )
    parser.add_argument(
        "--shell_path",
        required=False,
        help="Path of shell. Would be bash/dash/ash/sh",
        nargs="?",
        type=str,
        action="store",
        default=None,
    )
    parser.add_argument(
        "--name",
        required=False,
        help="Submission name",
        nargs="?",
        type=str,
        action="store",
        default=DEFAULT_SUBMISSION_NAME,
    )
    parser.add_argument(
        "--description",
        required=False,
        help="Submission description",
        nargs="?",
        type=str,
        action="store",
        default=DEFAULT_SUBMISSION_DESCRIPTION,
    )
    parser.add_argument(
        "--cpu",
        required=False,
        help="Number of CPU to be used",
        nargs="?",
        type=float,
        action="store",
        default=DEFAULT_SUBMISSION_CPU,
    )
    parser.add_argument(
        "--mem",
        required=False,
        help="Number of Memory to be used, can be 1024-based (e.g., 1024KiB) or 1000-based (100KB)",
        nargs="?",
        type=Union[str, float],
        action="store",
        default=DEFAULT_SUBMISSION_MEM,
    )
    parser.add_argument(
        "--cwd", required=False, help="Working directory", nargs="?", type=str, action="store", default=os.getcwd()
    )
    parser.add_argument(
        "--stdin",
        required=False,
        help="File used to be stdin, can be None",
        nargs="?",
        type=str,
        action="store",
        default=None,
    )
    parser.add_argument(
        "--stdout",
        required=False,
        help="File used to be stdout, can be None",
        nargs="?",
        type=str,
        action="store",
        default=None,
    )
    parser.add_argument(
        "--stderr",
        required=False,
        help="File used to be stderr, can be None",
        nargs="?",
        type=str,
        action="store",
        default=None,
    )
    parser.add_argument(
        "--no_preserve_env", required=False, help="Preserve current environment variables", nargs="?", action="store"
    )
    parser.add_argument(
        "--additional_env",
        required=False,
        help="Environments to add or replace. Specified in [NAME]=[VALUE] format",
        nargs="*",
        action="store",
    )
    parser.add_argument(
        "--depends", required=False, help="Depend on the finish of which submission id", nargs="*", action="store"
    )
    parser.add_argument("--tags", required=False, help="Tags used in filtering, etc.", nargs="*", action="store")
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    cl = YSJSCluster(conn=args.connection)
    _lh.info("YSJS %s Cluster %s -- %s", cl.config.schedule_method, cl.config.name, cl.config.description)
    env = {}
    if not args.no_preserve_env:
        env.update(os.environ)
    if args.additional_env is not None:
        for env_kv in args.additional_env:
            env_kv: str
            equal_pos = env_kv.find("=")
            if equal_pos == -1:
                exit(1)
            env[env_kv[:equal_pos]] = env_kv[equal_pos + 1 :]
    submission = YSJSSubmission.new(
        submission_name=args.name,
        submission_description=args.description,
        cpu=args.cpu,
        mem=1024 * 1024,  # TODO
        cwd=args.cwd,
        stdin=args.stdin,
        stdout=args.stdout,
        stderr=args.stderr,
        script_path=args.script_path,
        shell_path=args.shell_path,
        env=env,
        tags=args.tags,
        depends=args.depends,
    )
    job_id = cl.submit(submission)

    _lh.info(
        "Submission %s (id: %s, time: %s) success -> job_id %d",
        submission.submission_name,
        submission.submission_id,
        str(datetime.datetime.fromtimestamp(submission.submission_time)),
        job_id,
    )
