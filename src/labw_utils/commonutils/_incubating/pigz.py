import gzip
import multiprocessing
import operator
from functools import reduce
from typing import Callable, TypeVar, Iterable

import joblib
import tqdm

_InType = TypeVar("_InType")
_OutType = TypeVar("_OutType")


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
    :param backend: The backend to be used. Recommended to use ``loky``. You may also try ``threading`` if ``loky`` fails.
    :return: Generated new iterable.
    """
    it: Iterable[_OutType] = joblib.Parallel(
        n_jobs=n_jobs, backend=backend
    )(
        joblib.delayed(f)(i) for i in input_iterable
    )
    return it


def window(seq: bytes, n):
    """
    >>> list(window(b"111", 1))
    [b'1', b'1', b'1']
    >>> list(window(b"111", 2))
    [b'11', b'1']
    >>> list(window(b"111", 3))
    [b'111']
    >>> list(window(b"111", 4))
    [b'111']
    """
    cur_pos = 0
    while True:
        if cur_pos + n > len(seq):
            yield seq[cur_pos:]
            return
        elif cur_pos + n == len(seq):
            yield seq[cur_pos: cur_pos + n]
            return
        else:
            yield seq[cur_pos: cur_pos + n]
            cur_pos += n


def compress_iterator(data: Iterable[bytes], level: int) -> Iterable[bytes]:
    return parallel_map(
        lambda datum: gzip.compress(datum, compresslevel=level),
        tqdm.tqdm(list(data))
    )


def compress(data: bytes, level: int) -> bytes:
    return reduce(
        operator.add,
        compress_iterator(window(data, 1024 * 1024), level)
    )


if __name__ == "__main__":
    data_len = 1024 * 1024 * 4
    with open("/dev/random", "rb") as random:
        data = random.read(data_len)

    assert gzip.decompress(compress(data, 9)) == data
