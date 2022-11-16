from typing import Dict, Any, List, Iterable, TypeVar

_InType = TypeVar("_InType")


def iterable_translate(in_iterable: Iterable[_InType], trans_dict: Dict[_InType, _InType]) -> Iterable[_InType]:
    for old_item in in_iterable:
        if old_item in trans_dict.keys():
            yield trans_dict[old_item]
        else:
            yield old_item


def dict_translate(in_dict: Dict[_InType, Any], trans_dict: Dict[_InType, _InType]) -> Dict[_InType, Any]:
    """
    Dictionary Translator.

    This function will change the key of ``in_dict`` with the rules specified
    in ``trans_dict``.

    For example:

    >>> dict_translate({'A':1, 'B':2, 'C':3}, {'A':'a', 'B':'b'})
    {'a': 1, 'b': 2, 'C': 3}

    .. warning::
     The order of item will change!

    :param in_dict: The input dictionary.
    :param trans_dict: The translator.
    """
    return {k: v for k, v in zip(iterable_translate(in_dict.keys(), trans_dict), in_dict.values())}


def list_translate(in_list: List[_InType], trans_dict: Dict[_InType, _InType]) -> List[_InType]:
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
