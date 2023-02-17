"""
General-purposed helpers for Numpy NDArray and Torch Tensor.
"""

__all__ = (
    "scale_np_array",
    "scale_torch_array",
    "describe"
)

try:
    import torch
    from labw_utils.mlutils._ndarray_helper_with_torch import *
except ImportError:
    torch = None
    from labw_utils.mlutils._ndarray_helper_without_torch import *
