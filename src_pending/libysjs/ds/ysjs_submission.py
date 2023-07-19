from __future__ import annotations

import os
import shutil
import time
import uuid

from labw_utils.typing_importer import Union, Mapping, Optional, Any, List, Iterable

DEFAULT_SUBMISSION_NAME = "Unnamed"
DEFAULT_SUBMISSION_DESCRIPTION = "No description"
DEFAULT_SUBMISSION_MEM = 1024 * 1024
DEFAULT_SUBMISSION_CPU = 1


class YSJSSubmission:
    _submission_id: str
    _submission_name: str
    _submission_description: str
    _cpu: Union[int, float]
    _mem: Union[int, float]
    _submission_time: float
    _cwd: str
    _tags: List[str]
    _env: Mapping[str, str]
    _stdin: Optional[str]
    _stdout: Optional[str]
    _stderr: Optional[str]
    _script_path: str
    _shell_path: str
    _depends: List[str]

    def __init__(
            self,
            submission_id: str,
            submission_name: str,
            submission_description: str,
            cpu: Union[int, float],
            mem: Union[int, float],
            submission_time: float,
            cwd: str,
            env: Mapping[str, str],
            stdin: Optional[str],
            stdout: Optional[str],
            stderr: Optional[str],
            script_path: str,
            shell_path: str,
            tags: List[str],
            depends: List[str]
    ):
        self._submission_id = submission_id
        self._submission_name = submission_name
        self._cpu = cpu
        self._mem = mem
        self._submission_description = submission_description
        self._submission_time = submission_time
        self._cwd = cwd
        self._env = env
        self._stdin = stdin
        self._stdout = stdout
        self._stderr = stderr
        self._shell_path = shell_path
        self._script_path = script_path
        self._tags = tags
        self._depends = depends

    @classmethod
    def new(
            cls,
            script_path: str,
            cpu: Union[int, float] = DEFAULT_SUBMISSION_CPU,
            mem: Union[int, float] = DEFAULT_SUBMISSION_MEM,
            submission_name: str = DEFAULT_SUBMISSION_NAME,
            submission_description: str = DEFAULT_SUBMISSION_DESCRIPTION,
            cwd: Optional[str] = None,
            env: Optional[Mapping[str, str]] = None,
            stdin: Optional[str] = None,
            stdout: Optional[str] = None,
            stderr: Optional[str] = None,
            shell_path: Optional[str] = None,
            tags: Optional[List[str]] = None,
            depends: Optional[List[str]] = None
    ):
        if tags is None:
            tags = []
        if cwd is None:
            cwd = os.getcwd()
        if depends is None:
            depends = []
        if env is None:
            env = dict(os.environ)
        if shell_path is None:
            shell_path = shutil.which("bash")
        if shell_path is None:
            shell_path = shutil.which("dash")
        if shell_path is None:
            shell_path = shutil.which("ash")
        if shell_path is None:
            shell_path = shutil.which("sh")
        if shell_path is None:
            raise ValueError(
                "Cannot find suitable Shell; tried bash, dash, ash, sh"
            )
        script_path = os.path.abspath(script_path)
        cwd = os.path.abspath(cwd)
        if isinstance(stdin, str):
            stdin = os.path.abspath(stdin)
        if isinstance(stdout, str):
            stdout = os.path.abspath(stdout)
        if isinstance(stderr, str):
            stderr = os.path.abspath(stderr)

        return cls(
            submission_id=str(uuid.uuid4()),
            submission_name=submission_name,
            cpu=cpu,
            mem=mem,
            submission_description=submission_description,
            submission_time=time.time(),
            cwd=cwd,
            env=env,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            script_path=script_path,
            shell_path=shell_path,
            tags=tags,
            depends=depends
        )

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "submission_id": self._submission_id,
            "submission_name": self._submission_name,
            "cpu": self._cpu,
            "mem": self._mem,
            "submission_description": self._submission_description,
            "submission_time": self._submission_time,
            "cwd": self._cwd,
            "env": self._env,
            "stdin": self._stdin,
            "stdout": self._stdout,
            "stderr": self._stderr,
            "shell_path": self._shell_path,
            "script_path": self._script_path,
            "tags": self._tags,
            "depends": self._depends
        }

    def have_tag(self, tag: str) -> bool:
        return tag in self._tags

    @property
    def tags(self) -> Iterable[str]:
        return iter(self._tags)

    @property
    def submission_time(self) -> float:
        return self._submission_time

    @property
    def submission_id(self) -> str:
        return self._submission_id

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

    @property
    def cwd(self) -> str:
        return self._cwd

    @property
    def env(self) -> Mapping[str, str]:
        return self._env

    @property
    def stdin(self) -> Optional[str]:
        return self._stdin

    @property
    def stdout(self) -> Optional[str]:
        return self._stdout

    @property
    def stderr(self) -> Optional[str]:
        return self._stderr

    @property
    def script_path(self) -> str:
        return self._script_path

    @property
    def shell_path(self) -> str:
        return self._shell_path

    @property
    def cpu(self) -> Union[int, float]:
        return self._cpu

    @property
    def mem(self) -> Union[int, float]:
        return self._mem

    @property
    def submission_name(self) -> str:
        return self._submission_name

    @property
    def depends(self) -> List[str]:
        return self._depends
