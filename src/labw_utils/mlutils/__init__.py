"""
labw_utils.mlutils -- General-purposed machine- and deep-learning utilities.
"""

from labw_utils import UnmetDependenciesError
try:

    import numpy as np
    from numpy import typing as npt
except ImportError as e:
    raise UnmetDependenciesError("numpy") from e
