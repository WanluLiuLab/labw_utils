"""
typing_importer -- Python :py:mod:`typing`, :py:mod:`builtins` and :py:mod:`collections.abc` compatibility layer.

See also: PEP585

This also defines :py:class:`SequenceProxy` and :py:class:`MappingProxy` to create read-only view for underlying
:py:obj:`list` or :py:obj:`dict` objects.

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

An example on dicts:

>>> d = {1: 2, 3: 4, 5: 6}
>>> mv: MappingProxy[int, int] = MappingProxy(d)
>>> mv[1]
2
>>> len(mv)
3
>>> dict(mv)
{1: 2, 3: 4, 5: 6}
"""

from __future__ import annotations

__all__ = (
    "Callable",
    "Dict",
    "Tuple",
    "Set",
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

from typing import Any, Optional, Union, TypeVar, IO, TextIO, BinaryIO, AnyStr, NamedTuple, Generic

if sys.version_info >= (3, 9):
    from collections.abc import Callable, Iterable, Iterator, Awaitable, Coroutine, Generator, Mapping, \
        AsyncIterable, AsyncIterator, AsyncGenerator, Reversible, Container, Collection, MutableSet, MutableMapping, \
        Sequence, MutableSequence, ByteString, MappingView, KeysView, ItemsView, ValuesView, Sized
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
        KeysView, ItemsView, ValuesView, Sized

if sys.version_info < (3, 8):
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
    _seq: Sequence[_ItemType]

    def __getitem__(self, index: Union[int, slice]) -> _ItemType:
        return self._seq[index]

    def __len__(self) -> int:
        return len(self._seq)

    def __init__(self, seq: Iterable[_ItemType]):
        if isinstance(seq, Sequence):
            self._seq = seq
        elif isinstance(seq, Iterable):
            self._seq = list(seq)
        else:
            raise TypeError


class MappingProxy(Mapping[_KeyType, _ValueType]):
    _mapping: Mapping[_KeyType, _ValueType]

    def __getitem__(self, k: _KeyType) -> _ValueType:
        return self._mapping[k]

    def __len__(self) -> int:
        return len(self._mapping)

    def __iter__(self) -> Iterator[_KeyType]:
        yield from self._mapping

    def __init__(self, mapping: Mapping[_KeyType, _ValueType]):
        self._mapping = mapping
