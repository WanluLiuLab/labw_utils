"""
``labw_utils.bioutils.record.fai`` -- Record of Fasta Index (``*.fai``)

This data structure wraps one FAI record.

.. versionadded:: 1.0.2
"""

from __future__ import annotations

__all__ = ("FastaIndexRecordParserError", "FastaIndexRecord")

from labw_utils.devutils.decorators import create_class_init_doc_from_property


class FastaIndexRecordParserError(ValueError):
    """
    Error raised if:

    1. the FAI record have incorrect number of fields.
    2. the FAI record contains non-integer in integer fields.

    .. versionadded:: 1.0.2
    """

    def __init__(self, reason: str):
        super().__init__(reason)


@create_class_init_doc_from_property()
class FastaIndexRecord:
    """
    An entry from ``.fai`` files

    .. versionadded:: 1.0.2
    """

    _name: str
    _length: int
    _offset: int
    _line_blen: int
    _line_len: int

    __slots__ = ("_name", "_length", "_offset", "_line_blen", "_line_len")

    @property
    def name(self) -> str:
        """
        Chromosome Name
        """
        return self._name

    @property
    def length(self) -> int:
        """
        Chromosome length
        """
        return self._length

    @property
    def offset(self) -> int:
        """
        How far to `seek` to find this chromosome
        """
        return self._offset

    @property
    def line_blen(self) -> int:
        """
        Line length in bytes without newline
        """
        return self._line_blen

    @property
    def line_len(self) -> int:
        """
        Line length with newlines
        """
        return self._line_len

    def __init__(self, name: str, length: int, offset: int, line_len: int, line_blen: int):
        self._name = name
        self._length = length
        self._offset = offset
        self._line_blen = line_blen
        self._line_len = line_len

    @classmethod
    def from_fai_str(cls, fai_str: str):
        """
        Parse FAI string.

        :raises FastaIndexRecordParserError: On invalid input.
        """
        fields = fai_str.rstrip().split("\t")
        if len(fields) != 5:
            raise FastaIndexRecordParserError(f"Illegal record: '{fai_str}'. Need to have 5 fields.")
        try:
            new_instance = cls(
                name=fields[0],
                length=int(fields[1]),
                offset=int(fields[2]),
                line_blen=int(fields[3]),
                line_len=int(fields[4]),
            )
        except ValueError as e:
            raise FastaIndexRecordParserError(f"Illegal record: '{fai_str}'. Integer parsing error.") from e
        return new_instance

    def __repr__(self):
        return "\t".join((self.name, str(self.length), str(self.offset), str(self.line_blen), str(self.line_len)))

    def __str__(self):
        return repr(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FastaIndexRecord):
            raise TypeError(f"Object {other} (class: {type(other)}) is not FastaIndexRecord!")
        return repr(self) == repr(other)
