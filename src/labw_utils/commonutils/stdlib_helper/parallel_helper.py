"""
``labw_utils.stdlib_helper.parallel_helper`` -- Helper for Multiprocessing

This includes a very basic job pool and some helper classes

.. versionadded:: 1.0.2
"""

from __future__ import annotations

__all__ = (
    "PRIMITIVE_JOB_TYPE",
    "PROCESS_TYPE",
    "Job",
    "ParallelJobExecutor",
    "TimeOutKiller",
    "parallel_map",
    "easyexec",
)

import gc
import multiprocessing
import os
import subprocess
import threading
import time

from labw_utils import UnmetDependenciesError
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.devutils.decorators import create_class_init_doc_from_property
from labw_utils.typing_importer import Iterable
from labw_utils.typing_importer import Union, Optional, TypeVar, Callable, List

try:
    import joblib
except ImportError:
    raise UnmetDependenciesError("joblib")

PRIMITIVE_JOB_TYPE = Union[multiprocessing.Process, threading.Thread]
"""Python built-in job types"""

PROCESS_TYPE = None
"""
Process type that have a ``pid`` attribute.
Default is :py:class:`multiprocessing.Process`, :py:class:`subprocess.Popen`

If :external+psutil:py:mod:`psutil` is installed, will add :external+psutil:py:class:`psutil.Process`.
"""

try:
    import psutil

    PROCESS_TYPE = Union[multiprocessing.Process, subprocess.Popen, psutil.Process]
except ImportError:
    psutil = None
    PROCESS_TYPE = Union[multiprocessing.Process, subprocess.Popen]

_TERMINATE_HANDLER_TYPE = Callable[[PRIMITIVE_JOB_TYPE], None]
_CALLBACK_TYPE = Callable[[PRIMITIVE_JOB_TYPE], None]

_InType = TypeVar("_InType")
_OutType = TypeVar("_OutType")

_lh = get_logger(__name__)


@create_class_init_doc_from_property(
    text_after="""
:param terminate_handler: Function to terminate a job.
  Should accept the type of ``job_object`` and return :py:obj:`None`.
:param callback: Function executed after exit (successful or not) of a job.
  Should accept the type of ``job_object`` and return :py:obj:`None`.
"""
)
class Job:
    """
    The basic allocatable job with callback and termination function.
    It does not provide limitations on resource allocation or functions that handles errors.

    This class is NOT a subclass of :py:class:`threading.Thread`. It would not start new threads.
    All actions are synchronous.

    .. versionadded:: 1.0.2
    """

    _job_id: int
    _job_object: PRIMITIVE_JOB_TYPE
    _terminate_handler: Optional[_TERMINATE_HANDLER_TYPE]
    _callback: Optional[_CALLBACK_TYPE]

    def __init__(
        self,
        job_id: int,
        job_object: PRIMITIVE_JOB_TYPE,
        terminate_handler: Optional[_TERMINATE_HANDLER_TYPE] = None,
        callback: Optional[_CALLBACK_TYPE] = None,
    ):
        self._job_id = job_id
        self._job_object = job_object
        self._terminate_handler = terminate_handler
        self._callback = callback

    def start(self):
        """
        Start the job. This step should NOT raise exception.
        """
        self._job_object.start()

    def join(self):
        """
        Join the job and execute callback function if present. This step should NOT raise exception.
        """
        self._job_object.join()
        if self._callback is not None:
            self._callback(self._job_object)

    def terminate(self):
        """
        Terminate the job and execute callback, if present.
        """
        if self._terminate_handler is not None:
            self._terminate_handler(self._job_object)
        if self._callback is not None:
            self._callback(self._job_object)

    def terminate_without_callback(self):
        """
        Terminate the job without executing callback.
        """
        if self._terminate_handler is not None:
            self._terminate_handler(self._job_object)

    def __hash__(self):
        return self._job_id

    @property
    def job_id(self) -> int:
        """
        Unique job-id, which would also become ``hash`` of the job instance.
        """
        return self._job_id

    @property
    def job_object(self) -> PRIMITIVE_JOB_TYPE:
        """
        The internal object which a Job controls.
        """
        return self._job_object


@create_class_init_doc_from_property(
    text_after="""
:param delete_after_finish: Whether to delete job instance and perform garbage collection after finish
:param show_tqdm: Whether to show a progress bar.
"""
)
class ParallelJobExecutor(threading.Thread):
    """
    This is a parallel job executor,
    for jobs in a format of :py:class:`multiprocessing.Process` or :py:class:`threading.Thread`.

    This queue is designed for "batch" jobs.
    That is, the user should append all jobs before they start the executor.

    This executor is designed for non-stated jobs.
    That is, the executor will NOT save the state of any job.

    This executor would start a new thread.

    Length of the executor is defined by total number of jobs.

    .. versionadded:: 1.0.2
    """

    _pool_size: Union[int, float]
    _pool_name: str
    _refresh_interval: float

    _is_terminated: bool
    """
    Whether termination signal was sent to this thread
    """

    _is_appendable: bool
    """
    Whether this queue is appendable.
    """

    _pending_job_queue: List[Job]
    """
    Job waiting to be executed
    """

    _running_job_queue: List[Job]
    _finished_job_queue: List[Job]

    _n_jobs: int
    """
    Number of jobs to be executed
    """

    _delete_after_finish: bool
    _show_tqdm: bool

    def __init__(
        self,
        pool_name: str = "Unnamed pool",
        pool_size: Union[int, float] = 0,
        refresh_interval: float = 0.01,
        delete_after_finish: bool = True,
        show_tqdm: bool = True,
    ):
        super().__init__()
        self._pool_size = pool_size
        if self._pool_size == 0:
            self._pool_size = multiprocessing.cpu_count()
        self._pool_name = pool_name
        self._is_terminated = False
        self._is_appendable = True
        self._pending_job_queue = []
        self._running_job_queue = []
        self._finished_job_queue = []
        self._n_jobs = 0
        self._refresh_interval = refresh_interval
        self._delete_after_finish = delete_after_finish
        self._show_tqdm = show_tqdm

    def run(self):
        """
        Run the queue.
        """
        self._is_appendable = False
        if self._show_tqdm:
            pbar = tqdm(desc=self._pool_name, total=self._n_jobs)

        def _scan_through_process():
            """
            Scan through all processes and terminate the exited process.
            """
            for process in self._running_job_queue:
                if not process.job_object.is_alive():
                    process.join()
                    if isinstance(process.job_object, multiprocessing.Process):
                        process.job_object.close()
                    self._running_job_queue.remove(process)
                    if self._delete_after_finish:
                        del process
                        gc.collect()
                    else:
                        self._finished_job_queue.append(process)
                    if self._show_tqdm:
                        pbar.update(1)

        while len(self._pending_job_queue) > 0 and not self._is_terminated:
            while (
                len(self._pending_job_queue) > 0
                and len(self._running_job_queue) < self._pool_size
            ):
                new_process = self._pending_job_queue.pop(0)
                self._running_job_queue.append(new_process)
                new_process.start()
            _scan_through_process()
            time.sleep(self._refresh_interval)
        while len(self._running_job_queue) > 0 and not self._is_terminated:
            _scan_through_process()
            time.sleep(self._refresh_interval)
        self._is_terminated = True

    def stop(self):
        """
        Send termination signal.
        This will stop the job queue from adding more jobs.
        It would not affect job that is being executed.
        """
        self._is_appendable = False
        self._is_terminated = True
        self._pending_job_queue.clear()
        for job in self._running_job_queue:
            job.terminate()

    def append(
        self,
        mp_instance: PRIMITIVE_JOB_TYPE,
        terminate_handler: Optional[_TERMINATE_HANDLER_TYPE] = None,
        callback: Optional[_CALLBACK_TYPE] = None,
    ):
        """
        Commit a new job to the queue
        """
        if self._is_appendable:
            self._pending_job_queue.append(
                Job(
                    job_object=mp_instance,
                    job_id=self._n_jobs,
                    terminate_handler=terminate_handler,
                    callback=callback,
                )
            )
            self._n_jobs += 1
        else:
            raise ValueError("Job queue not appendable!")

    def iter_running_jobs(self) -> Iterable[PRIMITIVE_JOB_TYPE]:
        for job in self._running_job_queue:
            yield job.job_object

    def iter_pending_jobs(self) -> Iterable[PRIMITIVE_JOB_TYPE]:
        for job in self._pending_job_queue:
            yield job.job_object

    def iter_finished_jobs(self) -> Iterable[PRIMITIVE_JOB_TYPE]:
        for job in self._finished_job_queue:
            yield job.job_object

    @property
    def pool_size(self) -> float:
        """
        How many jobs is allowed to be executed in one time.

        Use :py:obj:`math.inf` to set unlimited or ``0`` to auto determine.
        """
        return self._pool_size

    @property
    def refresh_interval(self) -> float:
        """
        Interval for probing job status.
        """
        return self._refresh_interval

    @property
    def pool_name(self) -> str:
        """
        name of pool to be showed on progress bar, etc.
        """
        return self._pool_name

    @property
    def num_running_jobs(self) -> int:
        return len(self._running_job_queue)

    @property
    def num_pending_jobs(self) -> int:
        return len(self._pending_job_queue)

    @property
    def num_finished_jobs(self) -> int:
        return len(self._finished_job_queue)

    def __len__(self):
        return self._n_jobs


class TimeOutKiller(threading.Thread):
    """
    A timer that kills a process if time out is reached.



    After reaching the timeout, the killer will firstly send SIGTERM (15).
    If the process is alive after 3 seconds, it will send SIGKILL (9).

    .. warning :: This would not check whether the PIDs are in the same round.

    .. versionadded:: 1.0.2
    """

    _pid: int
    """
    Monitored process ID
    """

    _timeout: float

    def __init__(self, process_or_pid: Union[PROCESS_TYPE, int], timeout: float = 30.0):
        """
        .. warning :: Initialize the object after starting the monitored process!

        :param process_or_pid: A process represented either by some class with ``pid`` attribute or by its PID.
        :param timeout: The timeout in seconds, default 30.0.
        """
        super().__init__()
        if isinstance(process_or_pid, int):
            self._pid = process_or_pid
        else:
            self._pid = process_or_pid.pid
        self._timeout = timeout

    def run(self):
        time.sleep(self._timeout)
        try:
            os.kill(self._pid, 15)
        except (ProcessLookupError, PermissionError):
            return
        time.sleep(3)
        try:
            os.kill(self._pid, 9)
        except (ProcessLookupError, PermissionError):
            pass


def parallel_map(
    f: Callable[[_InType], _OutType],
    input_iterable: Iterable[_InType],
    n_jobs: int = multiprocessing.cpu_count(),
    backend: str = "threading",
) -> Iterable[_OutType]:
    """
    The parallel version of Python :py:func:`map` function (or, ``apply`` function in R)
    with :external+joblib:py:class:`joblib.Parallel` as backend.

    .. note ::
        This is to be refactored into a non :external+joblib:py:mod:`joblib`-dependent way.

    .. warning::
        With inappropriate parallelization, the system would consume lots of memory with minimal speed improvement!

    .. warning::
        Use with caution if you wish to parallely assign elements to an array.

    :param f: Function to be applied around an iterable.
    :param input_iterable: Iterable where a function would be applied to.
    :param n_jobs: Number of parallel threads. Would be max available CPU number if not set.
    :param backend: The backend to be used. Recommended to use ``threading``.
    :return: Generated new iterable.

    .. versionadded:: 1.0.2
    """
    it: Iterable[_OutType] = joblib.Parallel(n_jobs=n_jobs, backend=backend)(
        joblib.delayed(f)(i) for i in input_iterable
    )
    return it


def easyexec(
    cmd: List[str],
    log_path: Optional[str] = None,
    capture_output_path: Optional[str] = None,
    is_binary: bool = True,
    close_io: bool = True,
    raise_on_error: bool = True,
    cwd: str = os.getcwd(),
) -> int:
    _lh.debug(f"EASYEXEC {' '.join(cmd)} START")
    err_stream = (
        subprocess.DEVNULL if log_path is None else get_writer(log_path, is_binary)
    )
    out_stream = (
        err_stream
        if capture_output_path is None
        else get_writer(capture_output_path, is_binary)
    )
    p = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdin=subprocess.DEVNULL,
        stdout=out_stream,
        stderr=err_stream,
    )
    retv = p.wait()
    if close_io:
        try:
            out_stream.close()
        except AttributeError:
            pass
        try:
            err_stream.close()
        except AttributeError:
            pass
    if retv != 0:
        if raise_on_error:
            raise RuntimeError(f"EASYEXEC {' '.join(cmd)} FIN RETV={retv}")
        else:
            _lh.warning(f"EASYEXEC {' '.join(cmd)} FIN RETV={retv}")
    else:
        _lh.debug(f"EASYEXEC {' '.join(cmd)} FIN RETV={retv}")
    return retv
