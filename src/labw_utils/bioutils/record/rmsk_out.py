"""
``labw_utils.bioutils.record.rmsk_out`` -- Record of RepeatMasker Output files (``*.out``)

This data structure wraps one RMSK Out record.

.. versionadded:: 1.0.3
"""

from __future__ import annotations

__all__ = ("parse_record", "RmskOutRecord")

import re

from labw_utils.devutils.decorators import create_class_init_doc_from_property


def parse_record(record: str) -> RmskOutRecord:
    record_fields = re.sub(r" +", " ", record.strip()).split(" ")
    strand = record_fields[8] == "+"
    return RmskOutRecord(
        sw_score=int(record_fields[0]),
        pdiv=float(record_fields[1]),
        pdel=float(record_fields[2]),
        pins=float(record_fields[3]),
        qname=record_fields[4],
        qstart=int(record_fields[5].lstrip("(").rstrip(")")),
        qend=int(record_fields[6].lstrip("(").rstrip(")")),
        qleft=int(record_fields[7].lstrip("(").rstrip(")")),
        strand=strand,
        sname=record_fields[9],
        sfamily=record_fields[10],
        sstart=int(record_fields[11 if strand else 13].lstrip("(").rstrip(")")),
        send=int(record_fields[12].lstrip("(").rstrip(")")),
        sleft=int(record_fields[13 if strand else 11].lstrip("(").rstrip(")")),
        match_id=int(record_fields[14]),
        ssupressed=record.strip().endswith("*"),
    )


@create_class_init_doc_from_property()
class RmskOutRecord:
    """
    An entry from ``.fai`` files

    .. versionadded:: 1.0.2
    """

    #    SW   perc perc perc  query     position in query              matching       repeat                          position in repeat
    # score   div. del. ins.  sequence  begin    end          (left)   repeat         class/family                begin   end    (left)     ID
    #   307    0.0  0.0  0.0  chrX             1      262 (17718680) + (CTAAGC)n      Simple_repeat                     1    262     (0)     1

    _sw_score: int
    _pdiv: float
    _pdel: float
    _pins: float
    _qname: str
    _qstart: int
    _qend: int
    _qleft: int
    _sname: str
    _sfamily: str
    _sstart: int
    _send: int
    _sleft: int
    _ssupressed: bool
    _strand: bool
    _match_id: int

    __slots__ = (
        "_sw_score",
        "_pdiv",
        "_pdel",
        "_pins",
        "_qname",
        "_qstart",
        "_qend",
        "_qleft",
        "_sname",
        "_sfamily",
        "_sstart",
        "_send",
        "_sleft",
        "_strand",
        "_ssupressed",
        "_match_id",
    )

    @property
    def sw_score(self) -> int:
        """
        Smith-Waterman score of the match, usually complexity adjusted
        The SW scores are not always directly comparable.
        Sometimes the complexity adjustment has been turned off,
        and a variety of scoring-matrices are used dependent on repeat age and GC level.
        """
        return self._sw_score

    @property
    def pdiv(self) -> float:
        """% divergence = mismatches/(matches+mismatches)"""
        return self._pdiv

    @property
    def pdel(self) -> float:
        """% of bases opposite a gap in the query sequence (deleted bp)"""
        return self._pdel

    @property
    def pins(self) -> float:
        """% of bases opposite a gap in the repeat consensus (inserted bp)"""
        return self._pins

    @property
    def qname(self) -> str:
        """name of query sequence"""
        return self._qname

    @property
    def qstart(self) -> int:
        """starting position of match in query sequence"""
        return self._qstart

    @property
    def qend(self) -> int:
        """ending position of match in query sequence"""
        return self._qend

    @property
    def qleft(self) -> int:
        """no. of bases in query sequence past the ending position of match"""
        return self._qleft

    @property
    def sname(self) -> str:
        """name of the matching interspersed repeat"""
        return self._sname

    @property
    def sfamily(self) -> str:
        """the class of the repeat"""
        return self._sfamily

    @property
    def sleft(self) -> int:
        """no. of bases in repeat consensus sequence past the ending position of match"""
        return self._sleft

    @property
    def ssupressed(self) -> bool:
        """
        An asterisk (*) following the final column indicates that
        there is a higher-scoring match whose domain partly (<80%) includes the domain of the current match.
        """
        return self._ssupressed

    @property
    def sstart(self) -> int:
        """starting position of match in repeat consensus sequence"""
        return self._sstart

    @property
    def send(self) -> int:
        """ending position of match in repeat consensus sequence"""
        return self._send

    @property
    def strand(self) -> bool:
        """False: match is with the Complement of the repeat consensus sequence"""
        return self._strand

    @property
    def qlen(self) -> int:
        return self.qend + self.qleft

    @property
    def slen(self) -> int:
        return self.send + self.sleft

    @property
    def aln_dice(self) -> float:
        return (self.qend - self.qstart + self.send - self.sstart) / (self.qlen + self.slen)

    @property
    def paln(self) -> float:
        return (self.send - self.sstart) / self.slen

    @property
    def match_id(self) -> int:
        return self._match_id

    def __init__(
        self,
        sname: str,
        sfamily: str,
        sstart: int,
        send: int,
        sleft: int,
        ssupressed: bool,
        qname: str,
        qstart: int,
        qend: int,
        qleft: int,
        sw_score: int,
        pdiv: float,
        pdel: float,
        pins: float,
        strand: bool,
        match_id: int,
    ):
        self._sname = sname
        self._sfamily = sfamily
        self._sstart = sstart
        self._send = send
        self._sleft = sleft
        self._ssupressed = ssupressed
        self._qname = qname
        self._qstart = qstart
        self._qend = qend
        self._qleft = qleft
        self._sw_score = sw_score
        self._pdiv = pdiv
        self._pdel = pdel
        self._pins = pins
        self._strand = strand
        self._match_id = match_id
