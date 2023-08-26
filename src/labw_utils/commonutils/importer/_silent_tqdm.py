"""
A silent tqdm that does not pollutes stderr
"""

import sys

from labw_utils.devutils.decorators import create_class_init_doc_from_property
from labw_utils.typing_importer import Optional, Iterable, Sized, TypeVar, Iterator

__all__ = ("tqdm",)

_VarType = TypeVar("_VarType")


@create_class_init_doc_from_property()
class tqdm(Iterable[_VarType]):
    """
    A silent tqdm that does not pollute stderr.

    This module will be imported by :py:mod:`commonutils.importer.tqdm_importer` if stderr is not a TTY,
    which usually happens when the stderr is connected to a logger.

    This module does not support rich functions as TQDM.
    It only supports 4 quarters (i.e., it will change only 4 times)

    Examples:

    Using tqdm as an iterator:

    >>> list(tqdm(iterable=range(10)))
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    Using tqdm as context manager:

    >>> out_list = []
    >>> with tqdm(total=10) as pbar:
    ...    for i in range(10):
    ...        out_list.append(i)
    ...        pbar.update(1)
    >>> out_list
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    .. versionadded:: 1.0.2
    """

    __slots__ = ("_iterable", "_total", "_desc", "_n", "_quarters")

    _iterable: Optional[Iterable]
    _desc: Optional[str]
    _total: Optional[int]
    _n: int
    _quarters: float

    @property
    def total(self) -> Optional[int]:
        """
        The maximum number of items to be iterated.

        Set this variable if ``iterable`` is not sizeable or is ``None``.
        """
        return self._total

    @property
    def desc(self) -> Optional[str]:
        """
        Description show at the front of the progress bar.
        """
        return self._desc

    @property
    def iterable(self) -> Optional[Iterable[_VarType]]:
        """
        The iterable to be iterated when using ``tqdm`` in ``for`` loops.
        """
        return self._iterable

    def __iadd__(self, other: int):
        self.update(other)

    def __init__(
        self,
        iterable: Optional[Iterable[_VarType]] = None,
        desc: Optional[str] = None,
        total: Optional[int] = None,
        **kwargs,
    ):
        _ = kwargs  # Stop PyCharm from complaining
        self._desc = desc
        self._n = 0
        self._quarters = 0
        if total is None and iterable is not None:
            if isinstance(iterable, Sized):
                try:
                    self._total = len(iterable)
                except (TypeError, AttributeError):
                    self._total = None
            else:
                self._total = None
        else:
            self._total = total
        if total == float("inf"):
            self._total = None
        self._iterable = iterable

    def __iter__(self) -> Iterator[_VarType]:
        if self._iterable is None:
            raise TypeError
        for item in self._iterable:
            yield item
            self.update(1)

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, *args, **kwargs):
        return

    def update(self, i: int = 1):
        """
        Update the progress bar with given value.

        The given value should NOT be negative.
        """
        self._n += i
        if self._total:
            percent = round(self._n / self._total, 2)
            if percent > self._quarters:
                total_len = 100
                pbar_fill = "=" * int(self._quarters * total_len)
                pbar_blank = " " * int((1 - self._quarters) * total_len)
                print(
                    f"{self._desc}: {int(percent * 100)}% [{pbar_fill}|{pbar_blank}]",
                    file=sys.stderr,
                )
                self._quarters += 0.25
