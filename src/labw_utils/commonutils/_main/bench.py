"""
Benchmark Python interpreter.

.. warning::
    Not finished -- do not use

.. versionadded:: 1.0.2
"""
import concurrent.futures
import functools
import gc
import multiprocessing
import random
import time

try:
    import cProfile as profile
except ImportError:
    import profile


def recursive_fibonacci(n: int) -> int:
    def _recursive_fibonacci(_n: int) -> int:
        if _n == 1 or _n == 2:
            return 1
        else:
            return _recursive_fibonacci(_n - 1) + _recursive_fibonacci(_n - 2)

    return _recursive_fibonacci(n)


def recursive_fibonacci_cached(n: int) -> int:
    @functools.lru_cache(maxsize=4096)
    def _recursive_fibonacci(_n: int) -> int:
        if _n == 1 or _n == 2:
            return 1
        else:
            return _recursive_fibonacci(_n - 1) + _recursive_fibonacci(_n - 2)

    return _recursive_fibonacci(n)


def malloc_mfree_gc(n: int, size: int):
    for _ in range(n):
        m = bytearray(size)
        _ = m
        del m
        gc.collect()


def malloc_mfree(n: int, size: int):
    for _ in range(n):
        m = bytearray(size)
        _ = m
        del m


def _single_monte_caro_pi(_) -> bool:
    return random.uniform(-1, 1) ** 2 + random.uniform(-1, 1) ** 2 < 1


def monte_caro_pi(n: int) -> float:
    return sum(map(_single_monte_caro_pi, range(n))) / n * 4


def monte_caro_pi_threaded(n: int) -> float:
    with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as tpe:
        return sum(tpe.map(_single_monte_caro_pi, range(n))) / n * 4


def sqrt_binsearch(n: float, precision: float):
    lp, rp = 0.0, n
    pivot = (lp + rp) / 2
    while rp - lp > precision:
        if pivot**2 > n:
            rp = pivot
        else:
            lp = pivot
        pivot = (lp + rp) / 2
    return pivot


def sqrt_babylonian(n: float, precision: float):
    x = 1
    y = n
    while abs(x - y) > precision:
        x = (x + y) / 2
        y = n / x
    return x


def main(_):
    with profile.Profile() as pf:
        print("init")
        t1 = time.time_ns()
        print(recursive_fibonacci(25))
        t2 = time.time_ns()
        print("init1", t2 - t1)
        print(recursive_fibonacci_cached(25))
        t3 = time.time_ns()
        print("init2", t3 - t2)
        print(sqrt_binsearch(2, 0.00001))
        print(sqrt_babylonian(2, 0.00001))
        print(monte_caro_pi(100000))
        print(monte_caro_pi_threaded(100000))
        malloc_mfree_gc(20, 1024 * 1024 * 1024)
        malloc_mfree(20 * 1024 * 1024, 1)
        pf.dump_stats("1.pstat")
