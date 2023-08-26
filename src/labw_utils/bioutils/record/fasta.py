"""
``labw_utils.bioutils.record.fasta`` -- An In-Memory Fasta Record

.. versionadded:: 1.0.2
"""

__all__ = ("FastaRecord",)

from labw_utils.devutils.decorators import create_class_init_doc_from_property


@create_class_init_doc_from_property()
class FastaRecord:
    """
    A naive in-memory FASTQ record.

    .. versionadded:: 1.0.2
    """

    __slots__ = ("_seq_id", "_sequence")
    _seq_id: str
    _sequence: str

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

    def __init__(self, seq_id: str, sequence: str):
        self._seq_id = seq_id
        self._sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return "\n".join(
            (
                f">{self._seq_id}",
                self._sequence,
            )
        )

    def __str__(self):
        return repr(self)
