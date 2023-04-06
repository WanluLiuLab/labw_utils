"""
fasta.py -- An In-Memory Fasta Record
"""

__all__ = (
    "FastaRecord"
)


class FastaRecord:
    """
    A naive in-memory FASTQ record.
    """

    __slots__ = (
        '_seq_id',
        '_sequence'
    )
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
        return "\n".join((
            f">{self._seq_id}",
            self._sequence,
        ))

    def __str__(self):
        return repr(self)
