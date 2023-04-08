"""
labw_utils.mlutils.torch_layers -- Additional pyTorch layers.
"""


from labw_utils import UnmetDependenciesError

try:
    import torch
except ImportError as e:
    raise UnmetDependenciesError("torch") from e


from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.mlutils.ndarray_helper import describe


class Describe(torch.nn.Module):
    """
    The Describe Layer of PyTorch Module.

    Prints the description of matrix generated from last layer and pass the matrix without modification.
    """

    def __init__(self, prefix: str = ""):
        """
        The initializer

        :param prefix: Prefix of the printed message. Recommended to be the name of previous layer.

        See also: :py:func:`labw_utils.mlutils.ndarray_helper.describe`.
        """
        super().__init__()
        self._prefix = prefix
        self._lh = get_logger(__name__)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """"""
        self._lh.debug(self._prefix + describe(x))
        return x
