"""
.. versionadded:: 1.0.3
"""


import argparse
from abc import ABC, abstractmethod

__all__ = ("ArgparseRepresentableInterface",)


class ArgparseRepresentableInterface(ABC):
    """
    TODO
    """

    @staticmethod
    @abstractmethod
    def patch_agparser(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
        """Patch existing parser"""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_parsed_args(cls, parsed_args: argparse.Namespace):
        """Load the item from parsed arguments"""
        raise NotImplementedError
