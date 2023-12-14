"""
`labw_utils.stdlib_helper.itertools_helper` -- Some general-purposed functions on iterables.
    
.. versionadded:: 1.0.2
"""

from __future__ import annotations

__all__ = (
    "iterable_translate",
    "list_translate",
    "dict_translate",
    "window",
    "head",
    "tail",
)

from labw_utils.typing_importer import (
    Iterable,
    TypeVar,
    Mapping,
    List,
    Literal,
    Tuple,
    Optional,
)

_InType = TypeVar("_InType")
_VarType = TypeVar("_VarType")


def iterable_translate(in_iterable: Iterable[_InType], trans_dict: Mapping[_InType, _InType]) -> Iterable[_InType]:
    """
    Iterable translator.

    This function will change the elements of ``in_iterable`` with the rules specified in ``trans_dict``.

    .. seealso :: :py:func:`list_translate`.

    .. versionadded:: 1.0.2
    """
    trans_dict = dict(trans_dict)
    for old_item in in_iterable:
        if old_item in trans_dict.keys():
            yield trans_dict[old_item]
        else:
            yield old_item


def dict_translate(
    in_dict: Mapping[_InType, _VarType], trans_dict: Mapping[_InType, _InType]
) -> Mapping[_InType, _VarType]:
    """
    Dictionary Translator.

    This function will change the key of ``in_dict`` with the rules specified
    in ``trans_dict``.

    For example:

    >>> dict_translate({'A':1, 'B':2, 'C':3}, {'A':'a', 'B':'b'})
    {'a': 1, 'b': 2, 'C': 3}

    :param in_dict: The input dictionary.
    :param trans_dict: The translator.

    .. versionadded:: 1.0.2
    """
    trans_dict = dict(trans_dict)
    return {k: v for k, v in zip(iterable_translate(in_dict.keys(), trans_dict), in_dict.values())}


def list_translate(in_list: List[_InType], trans_dict: Mapping[_InType, _InType]) -> List[_InType]:
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

    .. versionadded:: 1.0.2
    """
    return list(iterable_translate(iter(in_list), trans_dict))


def window(
    it: Iterable[_InType],
    size: int,
    last_action: Literal["padd_front", "padd_back", "error", "ignore", "truncate"] = "ignore",
    padding: Optional[_InType] = None,
) -> Iterable[Tuple[_InType, ...]]:
    """
    Complex windowing support.

    >>> list(window([1, 2, 3], 2, "padd_front", 0))
    [(1, 2), (0, 3)]
    >>> list(window([1, 2, 3], 2, "padd_back", 0))
    [(1, 2), (3, 0)]
    >>> list(window([1, 2, 3], 2, "error", 0))
    Traceback (most recent call last):
        ...
    ValueError
    >>> list(window([1, 2, 3], 2, "ignore", 0))
    [(1, 2), (3,)]
    >>> list(window([1, 2, 3], 2, "truncate", 0))
    [(1, 2)]


    :param it: An iterable.
    :param size: Window size.
    :param last_action: What to do on insufficient last element.
    :param padding: What to pad if ``last_action`` is ``"padd_*"``.
    :return: Iterable of windows.

    .. versionadded:: 1.0.2
    """
    curr_list = []
    for item in it:
        curr_list.append(item)
        if len(curr_list) == size:
            yield tuple(curr_list)
            curr_list = []
    if len(curr_list) == size:
        yield tuple(curr_list)
    elif len(curr_list) == 0:
        return
    else:
        if last_action == "error":
            raise ValueError
        elif last_action == "truncate":
            return
        elif last_action == "padd_front":
            while len(curr_list) < size:
                curr_list.insert(0, padding)
        elif last_action == "padd_back":
            while len(curr_list) < size:
                curr_list.append(padding)
        yield tuple(curr_list)
        return


def head(it: Iterable[_VarType], n: int = 10) -> Iterable[_InType]:
    """
    Get the first ``n`` element of an iterable.

    .. versionadded:: 1.0.2
    """
    i = 0
    for item in it:
        i += 1
        if i == n:
            return
        yield item


def tail(it: Iterable[_VarType], n: int = 10) -> Iterable[_InType]:
    """
    Get the last ``n`` element of an iterable.

    .. versionadded:: 1.0.2
    """
    retl = []
    for item in it:
        if len(retl) < n:
            retl.append(item)
        else:
            retl.pop(0)
            retl.append(item)
    return iter(retl)
