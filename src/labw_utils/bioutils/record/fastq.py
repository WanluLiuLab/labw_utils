"""
``labw_utils.bioutils.record.fastq`` -- An In-Memory Fastq Record.

.. versionadded:: 1.0.2

Construct from 4 lines:

>>> fastq_seq = ["@A\\n", "AGCT\\n", "+\\n", "0000\\n"]
>>> fastq_record = FastqRecord.from_str(fastq_seq)
>>> fastq_record.seq_id
'A'
>>> fastq_record.sequence
'AGCT'
>>> len(fastq_record)
4
>>> fastq_record.quality
'0000'
>>> print(str(fastq_record))
@A
AGCT
+
0000

Construct from 1 string:

>>> fastq_record = FastqRecord.from_single_str('@A\\nAGCT\\n+\\n0000')
>>> print(str(fastq_record))
@A
AGCT
+
0000

The attributes in ``fastq_record`` are read-only. Example:

>>> fastq_record.seq_id = "5"
Traceback (most recent call last):
    ...
AttributeError: ...
"""

__all__ = (
    "FastqRecordParserError",
    "MisFormattedFastqRecordError",
    "SequenceQualityLengthMismatchError",
    "FastqRecord",
)

from labw_utils.devutils.decorators import create_class_init_doc_from_property, doc_del_attr
from labw_utils.typing_importer import List


@doc_del_attr(["__init__"])
class FastqRecordParserError(ValueError):
    """
    Generic Fastq parsing error.

    .. versionadded:: 1.0.2
    """

    pass


@doc_del_attr(["__init__"])
class MisFormattedFastqRecordError(FastqRecordParserError):
    """
    .. versionadded:: 1.0.2
    """

    def __init__(self, reason: str):
        super().__init__(reason)


@doc_del_attr(["__init__"])
class SequenceQualityLengthMismatchError(FastqRecordParserError):
    """
    Error raised when sequence length and quality mismatches.

    .. versionadded:: 1.0.2
    """

    def __init__(self, seq_id: str, sequence: str, quality: str):
        super().__init__(
            f"Illegal FASTQ record '{seq_id}': sequence '{sequence}' and quality '{quality}' length not equal."
        )


@create_class_init_doc_from_property(
    text_after="""

:raises SequenceQualityLengthMismatchError: On illegal record.
.. versionadded:: 1.0.2
"""
)
class FastqRecord:
    """
    A naive in-memory FASTQ record.
    """

    __slots__ = ("_seq_id", "_sequence", "_quality")
    _seq_id: str
    _sequence: str
    _quality: str

    @property
    def seq_id(self) -> str:
        """
        Sequence ID.
        """
        return self._seq_id

    @property
    def sequence(self) -> str:
        """
        The sequence.
        """
        return self._sequence

    @property
    def quality(self) -> str:
        """
        The corresponding quality, whose length should be equal to ``sequence``.
        """
        return self._quality

    def __init__(self, seq_id: str, sequence: str, quality: str, full_header: bool = True):
        if len(sequence) != len(quality):
            raise SequenceQualityLengthMismatchError(seq_id, sequence, quality)
        if full_header:
            self._seq_id = seq_id
        else:
            self._seq_id = seq_id.split(" ")[0].split("\t")[0]
        self._sequence = sequence
        self._quality = quality

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return "\n".join((f"@{self._seq_id}", self._sequence, "+", self._quality))

    def __str__(self):
        return repr(self)

    @classmethod
    def from_str(cls, lines: List[str], full_header: bool = True):
        """
        Generate from FASTQ sequence, 4 lines.

        This method is set to generate record from arbitrary :py:mod:`typing.TextIO` readers.

        :raises MisFormattedFastqRecordError: If the record is invalid.
        """
        if len(lines) != 4:
            raise MisFormattedFastqRecordError("Should get a 4-element aray representing 4 FASTQ lines.")
        l1, l2, l3, l4 = lines
        if not l1.startswith("@"):
            raise MisFormattedFastqRecordError(f"Line 1 {l1} should start with @")
        if not l3.startswith("+"):
            raise MisFormattedFastqRecordError(f"Line 3 {l3} should start with +")
        new_instance = cls(
            seq_id=l1[1:].rstrip("\n\r"),
            sequence=l2.rstrip("\n\r"),
            quality=l4.rstrip("\n\r"),
            full_header=full_header,
        )
        return new_instance

    @classmethod
    def from_single_str(cls, input_str: str):
        """
        Generate from FASTQ sequence, 1 line.

        :raises MisFormattedFastqRecordError: If the record is invalid.
        """
        return cls.from_str(input_str.splitlines())
