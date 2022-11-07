"""
Utilities that may be helpful.
"""

import math
from typing import Tuple

CoordinateType = Tuple[float, float]


def euclid_distance(d1: CoordinateType, d2: CoordinateType) -> float:
    """
    Get Euclid distance.

    >>> euclid_distance((0, 0), (3, 4))
    5.0
    """
    return math.sqrt((d1[0] - d2[0]) ** 2 + (d1[1] - d2[1]) ** 2)


def manhattan_distance(d1: CoordinateType, d2: CoordinateType) -> float:
    """
    Get manhattan distance.

    >>> manhattan_distance((0, 0), (3, 4))
    7
    """
    return abs(d1[0] - d2[0]) + abs(d1[1] - d2[1])
