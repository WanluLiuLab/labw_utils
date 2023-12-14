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

import collections.abc
import copy
import os
import typing

if os.getenv("SPHINX_BUILD") is not None:
    # Prevent Sphinx from building Python's internal docs
    __all__ = (
        "SequenceProxy",
        "MappingProxy",
    )
else:
    __all__ = (  # type: ignore
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
        "MappingProxy",
        "overload",
    )

import sys

Any = typing.Any
from typing import Optional as Optional

Union = typing.Union
TypeVar = typing.TypeVar
TextIO = typing.TextIO
IO = typing.IO
BinaryIO = typing.BinaryIO
AnyStr = typing.AnyStr
NamedTuple = typing.NamedTuple
Generic = typing.Generic
overload = typing.overload

if sys.version_info >= (3, 9):
    Callable = collections.abc.Callable
    Iterable = collections.abc.Iterable
    Iterator = collections.abc.Iterator
    Awaitable = collections.abc.Awaitable
    Coroutine = collections.abc.Coroutine
    Generator = collections.abc.Generator
    Mapping = collections.abc.Mapping
    AsyncIterable = collections.abc.AsyncIterable
    AsyncIterator = collections.abc.AsyncIterator
    AsyncGenerator = collections.abc.AsyncGenerator
    Reversible = collections.abc.Reversible
    Container = collections.abc.Container
    Collection = collections.abc.Collection
    MutableSet = collections.abc.MutableSet
    MutableMapping = collections.abc.MutableMapping
    MutableSequence = collections.abc.MutableSequence
    from collections.abc import Sequence as Sequence

    ByteString = collections.abc.ByteString
    MappingView = collections.abc.MappingView
    KeysView = collections.abc.KeysView
    ItemsView = collections.abc.ItemsView
    ValuesView = collections.abc.ValuesView
    Sized = collections.abc.Sized
    Hashable = collections.abc.Hashable

    Counter = collections.Counter
    OrderedDict = collections.OrderedDict
    ChainMap = collections.ChainMap
    List = list
    Dict = dict
    Set = set
    Tuple = tuple
    Type = type
    Deque = collections.deque
    DefaultDict = collections.defaultdict
else:
    Callable = typing.Callable
    Iterable = typing.Iterable
    Iterator = typing.Iterator
    Awaitable = typing.Awaitable
    Coroutine = typing.Coroutine
    Generator = typing.Generator
    Mapping = typing.Mapping
    AsyncIterable = typing.AsyncIterable
    AsyncIterator = typing.AsyncIterator
    AsyncGenerator = typing.AsyncGenerator
    Reversible = typing.Reversible
    Container = typing.Container
    Collection = typing.Collection
    MutableSet = typing.MutableSet
    MutableMapping = typing.MutableMapping
    MutableSequence = typing.MutableSequence
    ByteString = typing.ByteString
    MappingView = typing.MappingView
    KeysView = typing.KeysView
    ItemsView = typing.ItemsView
    ValuesView = typing.ValuesView
    Sized = typing.Sized
    Hashable = typing.Hashable

    Counter = typing.Counter
    OrderedDict = typing.OrderedDict
    ChainMap = typing.ChainMap
    List = typing.List
    from typing import Dict as Dict
    from typing import Sequence as Sequence

    Set = typing.Set
    Tuple = typing.Tuple
    Type = typing.Type
    Deque = typing.Deque
    DefaultDict = typing.DefaultDict
try:
    from typing import Literal as Literal
    from typing import Final as Final

except (AttributeError, ImportError):
    # For Python 3.7 only
    Final = typing._SpecialForm("Union", doc="")
    Literal = typing._SpecialForm("Union", doc="")

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

    def __getitem__(self, index):
        return self._seq[index]

    def __repr__(self):
        return repr(self._seq)

    def __str__(self):
        return repr(self)

    def __len__(self) -> int:
        return len(self._seq)

    def __init__(self, seq: Iterable[_ItemType], deep_copy: Optional[bool] = None):
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

    @classmethod
    def empty(cls):
        """
        Generate sequence proxy of empty list.
        """
        return cls([])


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

    def __init__(self, mapping: Mapping[_KeyType, _ValueType], deep_copy: Optional[bool] = None):
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

    def __repr__(self):
        return repr(self._mapping)

    def __str__(self):
        return repr(self)
