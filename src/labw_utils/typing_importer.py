"""
typing_importer -- Python :py:mod:`typing`, :py:mod:`builtins` and :py:mod:`collections.abc` compatibility layer.

See also: PEP585
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
    "Sized"
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
