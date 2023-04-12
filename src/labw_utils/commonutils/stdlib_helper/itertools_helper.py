"""
labw_utils.stdlib_helper.itertools_helper -- Some general-purposed functions on iterables.
"""

from __future__ import annotations

__all__ = (
    "iterable_translate",
    "list_translate",
    "dict_translate"
)

from typing import Any, TypeVar
from collections.abc import Iterable

_InType = TypeVar("_InType")
_VarType = TypeVar("_VarType")


def iterable_translate(in_iterable: Iterable[_InType], trans_dict: dict[_InType, _InType]) -> Iterable[_InType]:
    """
    Iterable translator.

    This function will change the elements of ``in_iterable`` with the rules specified in ``trans_dict``.

    See also :py:func:`list_translate`.
    """
    trans_dict = dict(trans_dict)
    for old_item in in_iterable:
        if old_item in trans_dict.keys():
            yield trans_dict[old_item]
        else:
            yield old_item


def dict_translate(in_dict: dict[_InType, _VarType], trans_dict: dict[_InType, _InType]) -> dict[_InType, _VarType]:
    """
    Dictionary Translator.

    This function will change the key of ``in_dict`` with the rules specified
    in ``trans_dict``.

    For example:

    >>> dict_translate({'A':1, 'B':2, 'C':3}, {'A':'a', 'B':'b'})
    {'a': 1, 'b': 2, 'C': 3}

    :param in_dict: The input dictionary.
    :param trans_dict: The translator.
    """
    trans_dict = dict(trans_dict)
    return {k: v for k, v in zip(iterable_translate(in_dict.keys(), trans_dict), in_dict.values())}


def list_translate(in_list: list[_InType], trans_dict: dict[_InType, _InType]) -> list[_InType]:
    """
    List Translator.

    Translate the list as is specified in py:func:`dict_translate`.

    The order of the item will NOT be changed.

    >>> list_translate(['A', 'B', 'C'], {'A':'a', 'B':'b'})
    ['a', 'b', 'C']

    :param in_list: Input list
    :param trans_dict: The translator.
    :type trans_dict: dict
    :return: Translated dictionary.
    """
    return list(iterable_translate(iter(in_list), trans_dict))
