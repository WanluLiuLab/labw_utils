"""
``labw_utils.mlutils.ndarray_helper`` -- General-purposed helpers for Numpy NDArray and Torch Tensor.

.. versionadded:: 1.0.2
"""

from __future__ import annotations

from labw_utils import UnmetDependenciesError

__all__ = ("scale_np_array", "scale_torch_array", "describe")

import numpy as np
import numpy.typing as npt

from labw_utils.typing_importer import Any, Union, Optional

try:
    import torch
except ImportError:
    torch = None

if torch is None:
    _Tensor = npt.NDArray
else:
    _Tensor = Union[npt.NDArray, torch.Tensor]


def _scale_impl(
    x: _Tensor,
    out_range: tuple[Union[int, float], Union[int, float]],
    domain: tuple[Union[int, float], Union[int, float]],
) -> _Tensor:
    if domain[1] == domain[0]:
        return x
    y = (x - (domain[1] + domain[0]) / 2) / (domain[1] - domain[0])
    return y * (out_range[1] - out_range[0]) + (out_range[1] + out_range[0]) / 2


def scale_np_array(
    x: npt.NDArray,
    domain: Optional[tuple[Union[int, float], Union[int, float]]] = None,
    out_range: tuple[Union[int, float], Union[int, float]] = (0, 1),
) -> npt.NDArray:
    """
    Scale a Numpy array to specific range.

    .. seealso :: :py:func:`scale_torch_array`

    Example:

    >>> scale_np_array(np.array([1,2,3,4,5]), out_range=(0, 1))
    array([0.  , 0.25, 0.5 , 0.75, 1.  ])

    .. versionadded:: 1.0.2
    """
    if domain is None:
        domain = np.min(x), np.max(x)
    return _scale_impl(x, out_range, domain)


if torch is not None:

    def scale_torch_array(
        x: torch.Tensor,
        domain: Optional[tuple[Union[int, float], Union[int, float]]] = None,
        out_range: tuple[Union[int, float], Union[int, float]] = (0, 1),
    ) -> torch.Tensor:
        """
        Scale a Torch array to specific range.

        .. seealso :: :py:func:`scale_np_array`

        Example:

        >>> scale_torch_array(torch.tensor(np.array([1,2,3,4,5])), out_range=(0, 1))
        tensor([0.0000, 0.2500, 0.5000, 0.7500, 1.0000])

        .. versionadded:: 1.0.2
        """
        if domain is None:
            domain = torch.min(x).item(), torch.max(x).item()
        return _scale_impl(x, out_range, domain)

else:

    def scale_torch_array(
        x: Any, domain: Optional[Any] = None, out_range: tuple[Union[int, float], Union[int, float]] = (0, 1)
    ) -> "torch.Tensor":
        _ = x, domain, out_range
        del x, domain, out_range
        raise UnmetDependenciesError("pytorch")


def describe(array: _Tensor) -> str:
    """
    Describe the array by data type, shape, quantiles or unique values.

    Example:

    >>> import sys, pytest
    >>> if sys.platform.startswith('win'):
    ...     pytest.skip('this doctest does not work on Windows')
    ...

    >>> describe(np.array([0, 0, 1, 1]))
    'ndarray[int64] with shape=(4,); uniques=[0 1]; mean=0.5'

    >>> describe(np.array(np.random.uniform(0, 21, size=12000), dtype=int))
    "ndarray[int64] with shape=(12000,); quantiles=['0.00', '5.00', '10.00', '15.00', '20.00']; mean=..."

    :param array: The Numpy array or pyTorch Tensor to be described.
    :return: Description of the array.

    .. versionadded:: 1.0.2
    """
    q = [0, 0.25, 0.5, 0.75, 1]
    _shape = tuple(array.shape)
    if torch is not None and isinstance(array, torch.Tensor):
        _unique = array.unique()
    else:
        _unique = np.unique(array)
    if len(_unique) > 10:
        try:
            if torch is not None and isinstance(array, torch.Tensor):
                array = array.float()
                _quantiles = list(map(lambda _q: f"{array.quantile(q=_q):.2f}", q))
            else:
                _quantiles = list(map(lambda _q: f"{_q:.2f}", np.quantile(array, q=q)))
        except (IndexError, RuntimeError) as e:
            _quantiles = [f"ERROR {e}"]
        _quantiles_str = f"quantiles={_quantiles}"
    else:
        _quantiles_str = f"uniques={_unique}"
    _mean = array.mean()
    return f"{type(array).__name__}[{array.dtype}] with shape={_shape}; {_quantiles_str}; mean={_mean}"


class DimensionMismatchException(ValueError):
    """
    A friendly representation of dimension mismatch.

    .. versionadded:: 1.0.2
    """

    def __init__(self, arr1: _Tensor, arr2: _Tensor, arr1_name: str = "arr1", arr2_name: str = "arr2"):
        super().__init__(
            f"Array {arr1_name} and {arr2_name}  dimension mismatch!\n"
            f"\twhere {arr1_name} is {describe(arr1)}\n"
            f"\twhere {arr2_name} is {describe(arr2)}\n"
        )
