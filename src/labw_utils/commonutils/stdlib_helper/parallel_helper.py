"""
parallel_helper.py -- Helper for Multiprocessing

This includes a very basic job pool and some helper classes
"""

import gc
import multiprocessing
import os
import subprocess
import threading
import time
from typing import Union, List, Optional, Callable, TypeVar, Iterable

import joblib

from labw_utils.commonutils.importer.tqdm_importer import tqdm

_JOB_TYPE = Union[multiprocessing.Process, threading.Thread]
_PROCESS_TYPE = Union[multiprocessing.Process, subprocess.Popen]
_TERMINATE_HANDLER_TYPE = Callable[[_JOB_TYPE], None]
_CALLBACK_TYPE = Callable[[_JOB_TYPE], None]

_InType = TypeVar("_InType")
_OutType = TypeVar("_OutType")


class Job:
    _job_id: int
    _job_object: _JOB_TYPE
    _terminate_handler: Optional[_TERMINATE_HANDLER_TYPE]
    _callback: Optional[_CALLBACK_TYPE]

    def __init__(
            self,
            job_id: int,
            job_object: _JOB_TYPE,
            terminate_handler: Optional[_TERMINATE_HANDLER_TYPE] = None,
            callback: Optional[_CALLBACK_TYPE] = None
    ):
        self._job_id = job_id
        self._job_object = job_object
        self._terminate_handler = terminate_handler
        self._callback = callback

    def start(self):
        self._job_object.start()

    def join(self):
        self._job_object.join()
        if self._callback is not None:
            self._callback(self._job_object)

    def terminate(self):
        if self._terminate_handler is not None:
            self._terminate_handler(self._job_object)
        if self._callback is not None:
            self._callback(self._job_object)

    def terminate_without_callback(self):
        if self._terminate_handler is not None:
            self._terminate_handler(self._job_object)

    def __hash__(self):
        return self._job_id

    @property
    def job_id(self) -> int:
        return self._job_id

    @property
    def job_object(self) -> _JOB_TYPE:
        return self._job_object


class ParallelJobExecutor(threading.Thread):
    """
    This is a parallel job executor,
    for jobs in a format of :py:class:`multiprocessing.Process` or :py:class:`threading.Thread`.

    This queue is designed for "batch" jobs.
    That is, the user should append all jobs before they start the executor.

    This executor is designed for non-stated jobs.
    That is, the executor will NOT save the state of any job.
    """

    _pool_size: Union[int, float]
    """
    How many jobs is allowed to be executed in one time.

    Use ``math.inf`` to set unlimited or ``0`` to auto determine.
    """

    _pool_name: str
    """
    name of pool to be showed on progress bar, etc.
    """

    _refresh_interval: float
    """
    Interval for probing job status.
    """

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
    """
    Whether to delete job instance after finish
    """

    _show_tqdm: bool

    def __init__(
            self,
            pool_name: str = "Unnamed pool",
            pool_size: Union[int, float] = 0,
            refresh_interval: float = 0.01,
            delete_after_finish: bool = True,
            show_tqdm: bool = True
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
            while len(self._pending_job_queue) > 0 and len(self._running_job_queue) < self._pool_size:
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
        """
        self._is_appendable = False
        self._is_terminated = True
        self._pending_job_queue.clear()
        for job in self._running_job_queue:
            job.terminate()

    def append(
            self,
            mp_instance: _JOB_TYPE,
            terminate_handler: Optional[_TERMINATE_HANDLER_TYPE] = None,
            callback: Optional[_CALLBACK_TYPE] = None
    ):
        """
        Commit a new job to the queue
        """
        if self._is_appendable:
            self._pending_job_queue.append(Job(
                job_object=mp_instance,
                job_id=self._n_jobs,
                terminate_handler=terminate_handler,
                callback=callback
            ))
            self._n_jobs += 1
        else:
            raise ValueError("Job queue not appendable!")

    def iter_running_jobs(self) -> Iterable[_JOB_TYPE]:
        for job in self._running_job_queue:
            yield job.job_object

    def iter_pending_jobs(self) -> Iterable[_JOB_TYPE]:
        for job in self._running_job_queue:
            yield job.job_object

    def iter_finished_jobs(self) -> Iterable[_JOB_TYPE]:
        for job in self._finished_job_queue:
            yield job.job_object

    @property
    def num_running_jobs(self) -> int:
        return len(self._running_job_queue)

    @property
    def num_pending_jobs(self) -> int:
        return len(self._pending_job_queue)

    @property
    def num_finished_jobs(self) -> int:
        return len(self._finished_job_queue)


class TimeOutKiller(threading.Thread):
    """
    A timer that kills a process if time out is reached.

    A process can be either represented using :py:class:`multiprocessing.Process` or :py:class`subprocess.Popen`,
    or by its PID.

    After reaching the timeout, the killer will firstly send SIGTERM (15).
    If the process is alive after 3 seconds, it will send SIGKILL (9).

    TODO: Check whether the PIDs are in the same round.
    """

    _pid: int
    """
    Monitored process ID
    """

    timeout: float
    """
    The timeout in seconds, default 30.0
    """

    def __init__(self, process_or_pid: Union[_PROCESS_TYPE, int], timeout: float = 30.0):
        """
        .. warning :: Initialize the object after starting the monitored process!
        """
        super().__init__()
        if isinstance(process_or_pid, int):
            self._pid = process_or_pid
        else:
            self._pid = process_or_pid.pid
        self.timeout = timeout

    def run(self):
        time.sleep(self.timeout)
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
    The parallel version of Python :external:py:func:`map` function (or, ``apply`` function in R).

    See also: :external+joblib:py:class:`joblib.Parallel`.

    .. warning::
        With inappropriate parallelization, the system would consume lots of memory with minimal speed improvement!

    .. warning::
        Use with caution if you wish to parallely assign elements to an array.

    :param f: Function to be applied around an iterable.
    :param input_iterable: Iterable where a function would be applied to.
    :param n_jobs: Number of parallel threads. Would be max available CPU number if not set.
    :param backend: The backend to be used. Recommended to use ``threading``.
    :return: Generated new iterable.
    """
    it: Iterable[_OutType] = joblib.Parallel(
        n_jobs=n_jobs, backend=backend
    )(
        joblib.delayed(f)(i) for i in input_iterable
    )
    return it
