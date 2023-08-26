"""
``labw_utils.bioutils.algorithm.utils`` -- Utilities that may be helpful.

 .. versionadded:: 1.0.2
"""
from __future__ import annotations

import math

from labw_utils.typing_importer import Tuple, List

CoordinateType = Tuple[float, float]


def euclid_distance(d1: CoordinateType, d2: CoordinateType) -> float:
    """
    Get Euclid distance.

    >>> euclid_distance((0, 0), (3, 4))
    5.0

     .. versionadded:: 1.0.2
    """
    return math.sqrt((d1[0] - d2[0]) ** 2 + (d1[1] - d2[1]) ** 2)


def manhattan_distance(d1: CoordinateType, d2: CoordinateType) -> float:
    """
    Get manhattan distance.

    >>> manhattan_distance((0, 0), (3, 4))
    7

     .. versionadded:: 1.0.2
    """
    return abs(d1[0] - d2[0]) + abs(d1[1] - d2[1])


def merge_intervals(arr: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    See: <https://www.geeksforgeeks.org/merging-intervals/>

     .. versionadded:: 1.0.2
    """
    # Sorting based on the increasing order
    # of the start intervals
    if not arr:
        return []
    arr = list(list(it) for it in arr)
    arr.sort(key=lambda x: x[0])

    # Stores index of last element
    # in output array (modified arr[])
    index = 0

    # Traverse all input Intervals starting from
    # second interval
    for i in range(1, len(arr)):
        # If this is not first Interval and overlaps
        # with the previous one, Merge previous and
        # current Intervals
        if arr[index][1] >= arr[i][0]:
            arr[index][1] = max(arr[index][1], arr[i][1])
        else:
            index = index + 1
            arr[index] = arr[i]
    retl = []
    for i in range(index + 1):
        retl.append(arr[i])
    return retl
