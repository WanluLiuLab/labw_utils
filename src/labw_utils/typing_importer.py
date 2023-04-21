"""
typing_importer -- Python :py:mod:`typing`, :py:mod:`builtins` and :py:mod:`collections.abc` compatibility layer.

This module serves as a compatibility layer between Python <= 3.8 and >= 3.9.
It would import from :py:mod:`typing` or :py:mod:`collections.abc` for generalized types
according to :pep:`585`.

This also defines :py:class:`SequenceProxy` and :py:class:`MappingProxy` to create read-only view for underlying
:py:obj:`list` or :py:obj:`dict` objects.

.. versionadded:: 1.0.1
"""

from __future__ import annotations
import copy

__all__ = (
    "Callable",
    "Dict",
    "Tuple",
    "Set",
    "Hashable",
    "Type",
    "Deque",
    "Iterable",
    "Iterator",
    "List",
    "Awaitable",
    "Coroutine",
    "Generator",
    "Mapping",
    "AsyncIterable",
    "AsyncIterator",
    "AsyncGenerator",
    "Reversible",
    "Container",
    "Collection",
    "MutableSet",
    "MutableMapping",
    "Sequence",
    "MutableSequence",
    "DefaultDict",
    "OrderedDict",
    "Counter",
    "ChainMap",
    "ByteString",
    "MappingView",
    "ItemsView",
    "KeysView",
    "ValuesView",
    "Optional",
    "Any",
    "Union",
    "Literal",
    "TypeVar",
    "IO",
    "TextIO",
    "BinaryIO",
    "AnyStr",
    "NamedTuple",
    "Final",
    "Generic",
    "Sized",
    "SequenceProxy",
    "MappingProxy"
)

import sys
from abc import abstractmethod

from typing import Any, Optional, Union, TypeVar, IO, TextIO, BinaryIO, AnyStr, NamedTuple, Generic, overload

if sys.version_info >= (3, 9):
    from collections.abc import Callable, Iterable, Iterator, Awaitable, Coroutine, Generator, Mapping, \
        AsyncIterable, AsyncIterator, AsyncGenerator, Reversible, Container, Collection, MutableSet, MutableMapping, \
        Sequence, MutableSequence, ByteString, MappingView, KeysView, ItemsView, ValuesView, Sized, Hashable
    from collections import Counter, OrderedDict, defaultdict, ChainMap, deque

    List = list
    Dict = dict
    Set = set
    Tuple = tuple
    Type = type
    Deque = deque
    DefaultDict = defaultdict
else:
    from typing import List, Dict, Set, Tuple, Type, Callable, Deque, Iterable, Iterator, Awaitable, Coroutine, \
        Generator, Mapping, AsyncIterable, AsyncIterator, AsyncGenerator, Reversible, Container, Collection, MutableSet, \
        MutableMapping, Sequence, MutableSequence, Counter, OrderedDict, DefaultDict, ChainMap, ByteString, MappingView, \
        KeysView, ItemsView, ValuesView, Sized, Hashable

if sys.version_info < (3, 8):
    # For Python 3.7 only
    class _Subscriptable:
        def __getitem__(self, item):
            ...


    Final = _Subscriptable()
    Literal = _Subscriptable()

else:
    from typing import Literal, Final

_ItemType = TypeVar("_ItemType")
_KeyType = TypeVar("_KeyType")
_ValueType = TypeVar("_ValueType")


class SequenceProxy(Sequence[_ItemType]):
    """
    A read-only proxy to iterables and sequences.

    Valid input types are:

    - :py:obj:`list`.
    - :py:obj:`typing.Sequence`.
    - :py:obj:`typing.Iterable`.

    >>> l = [1, 2, 3]
    >>> lv = SequenceProxy(l)
    >>> lv[1]
    2
    >>> len(lv)
    3
    >>> list(reversed(lv))
    [3, 2, 1]

    Works also on iterables:

    >>> li = iter(l)
    >>> lv: SequenceProxy[int] = SequenceProxy(li)
    >>> lv[1]
    2
    >>> len(lv)
    3
    >>> list(reversed(lv))
    [3, 2, 1]

    .. versionadded:: 1.0.1
    """

    _seq: Sequence[_ItemType]

    def __instancecheck__(self, instance: object):
        return isinstance(instance, Sequence)

    def __subclasscheck__(self, subclass):
        return issubclass(subclass, Sequence)

    @overload
    def __getitem__(self, index: int) -> _ItemType:
        return self._seq[index]

    @overload
    def __getitem__(self, index: slice) -> Sequence[_ItemType]:
        return self._seq[index]

    def __getitem__(self, index: Union[int, slice]) -> Union[Sequence[_ItemType], _ItemType]:
        return self._seq[index]

    def __len__(self) -> int:
        return len(self._seq)

    def __init__(
            self,
            seq: Iterable[_ItemType],
            deep_copy: Optional[bool] = None
    ):
        """
        Initialize the class with input sequence.

        :param seq: Input sequence.
            If sequence is :py:obj:`typing.Iterable`, would be persisted as :py:obj:`list`.

        :param deep_copy: Whether to copy the sequence into local cache.

            - :py:obj:`None`: Not perform copying. Proxy the sequence as-is.
            - :py:obj:`False`: Persist the input sequence as :py:obj:`list` and perform shallow copy.
            - :py:obj:`True`: Persist the input sequence as :py:obj:`list` and perform deep copy.
        :raises TypeError: If the input type is not :py:obj:`typing.Iterable` .
        """
        if isinstance(seq, Iterable):
            seq = list(seq)  # Force persistance over Iterable
        elif not isinstance(seq, Sequence):
            raise TypeError
        if deep_copy is None:
            self._seq = seq
        elif deep_copy is False:
            self._seq = copy.copy(list(seq))
        else:
            self._seq = copy.deepcopy(list(seq))


class MappingProxy(Mapping[_KeyType, _ValueType]):
    """
    A read-only proxy to iterables and sequences.

    Valid input types are:

    - :py:obj:`dict`.
    - :py:obj:`typing.Mapping`

    >>> d = {1: 2, 3: 4, 5: 6}
    >>> mv: MappingProxy[int, int] = MappingProxy(d)
    >>> mv[1]
    2
    >>> len(mv)
    3
    >>> dict(mv)
    {1: 2, 3: 4, 5: 6}

    .. versionadded:: 1.0.1
    """
    _mapping: Mapping[_KeyType, _ValueType]

    def __getitem__(self, k: _KeyType) -> _ValueType:
        return self._mapping[k]

    def __len__(self) -> int:
        return len(self._mapping)

    def __iter__(self) -> Iterator[_KeyType]:
        yield from self._mapping

    def __init__(
            self,
            mapping: Mapping[_KeyType, _ValueType],
            deep_copy: Optional[bool] = None
    ):
        """
        Load the class with input mapping.

        :param mapping: Input mapping.
        :param deep_copy: Whether to copy the mapping into local cache.

            - :py:obj:`None`: Not perform copying. Proxy the mapping as-is.
            - :py:obj:`False`: Persist the input mapping as :py:obj:`dict` and perform shallow copy.
            - :py:obj:`True`: Persist the input mapping as :py:obj:`dict` and perform deep copy.
        :raises TypeError: If the input type is not :py:obj:`typing.Mapping`.
        """
        if not isinstance(mapping, Mapping):
            raise TypeError
        if deep_copy is None:
            self._mapping = mapping
        elif deep_copy is False:
            self._mapping = copy.copy(dict(mapping))
        else:
            self._mapping = copy.deepcopy(dict(mapping))
