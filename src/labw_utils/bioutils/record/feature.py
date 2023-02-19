"""
feature.py -- General-Purposed GTF/GFF3/BED Record that Represents a Genomic Feature

This module includes GTF/GFF3/BED record datastructure and their one-line parsers.
"""

from __future__ import annotations

import enum
from abc import abstractmethod
from functools import total_ordering
from typing import Union, Optional, Dict, Iterable, List

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)

GtfAttributeValueType = Union[
    str, int, float, bool, None,
    List[str], List[int], List[float], List[bool], List[None]
]

GTFAttributeType = Dict[str, GtfAttributeValueType]
"""Type of GTF/GFF fields"""

VALID_GTF_QUOTE_OPTIONS = (
    "none",
    "blank",
    "string",
    "all"
)
"""
Valid GTF Quoting Options. They are:

* ``none`` Never quote, even if blanks found inside.
* ``blank``: Quote if blanks (\\r, \\n, \\t, space, etc.) found inside.
* ``string``: Quote if the field have type string. Will not quote for numeric types.
* ``all``: Quote all fields.
"""

DEFAULT_GTF_QUOTE_OPTIONS = "all"


class _NotSet:
    pass


_notset = _NotSet()


def feature_repr(v: GtfAttributeValueType) -> str:
    """
    Python standard :py:func:`repr` for genomic data.

    >>> feature_repr(None)
    '.'
    >>> feature_repr(1)
    '1'
    >>> feature_repr(1.0)
    '1.0'
    >>> feature_repr("ANS")
    'ANS'
    >>> feature_repr([])
    Traceback (most recent call last):
        ...
    TypeError: <class 'list'> is not supported!
    """
    if v is None:
        attr_str = "."
    elif isinstance(v, str):
        attr_str = v
    elif isinstance(v, float) or isinstance(v, int) or isinstance(v, bool):
        attr_str = str(v)
    else:
        raise TypeError(f"{type(v)} is not supported!")
    return attr_str


class FeatureParserError(ValueError):
    pass


class RegionError(FeatureParserError):
    def __init__(self, *args):
        super(RegionError, self).__init__(*args)


class FeatureType(enum.IntEnum):
    NOT_PRESENT = -1
    UNKNOWN = 0
    EXON = 9
    FIVE_PRIME_UTR = 8
    THREE_PRIME_UTR = 7
    OTHER_UTR = 6
    CDS = 5
    START_CODON = 4
    STOP_CODON = 3
    TRANSCRIPT = 2
    GENE = 1


_raw_feature_type_translator = {
    "3utr": FeatureType.THREE_PRIME_UTR,
    "three_prime_utr": FeatureType.THREE_PRIME_UTR,
    "5utr": FeatureType.FIVE_PRIME_UTR,
    "five_prime_utr": FeatureType.FIVE_PRIME_UTR,
    "utr": FeatureType.OTHER_UTR,
    "transcript": FeatureType.TRANSCRIPT,
    "gene": FeatureType.GENE,
    "exon": FeatureType.EXON,
    "cds": FeatureType.CDS,
    "start_codon": FeatureType.START_CODON,
    "stop_codon": FeatureType.STOP_CODON,
}


class BiologicalIntervalInterface:

    @property
    @abstractmethod
    def seqname(self) -> str:
        """
        Chromosome or Contig name.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def start(self) -> int:
        """
        Inclusive 1-based start position.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def start0b(self) -> int:
        """
        Inclusive 0-based start position.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def end(self) -> int:
        """
        Inclusive 1-based end position.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def end0b(self) -> int:
        """
        Exclusive 0-based end position.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def strand(self) -> Optional[bool]:
        """
        True (``+``) or False (``-``) or None (``.``)
        """
        raise NotImplementedError

    @abstractmethod
    def regional_equiv(self, other: BiologicalIntervalInterface):
        raise NotImplementedError


class FeatureInterface(BiologicalIntervalInterface):

    @property
    @abstractmethod
    def source(self) -> Optional[str]:
        """
        The source of this record. e.g. ``hg38_rmsk`` or ``ensembl``.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def feature(self) -> Optional[str]:
        """
        Feature type name. e.g. ``exon`` or ``start_codon`` or ``5UTR``.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def parsed_feature(self) -> FeatureType:
        raise NotImplementedError

    @property
    @abstractmethod
    def score(self) -> Optional[Union[int, float]]:
        """
        Some kind of scoring.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def frame(self) -> Optional[int]:
        """
        One of ``0`` (first base of the feature is the first base of a codon),
                        ``1`` (the second base is the first base of a codon) or ``2``.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def attribute_keys(self) -> Iterable[str]:
        """Other attributes presented in Key-Value pair"""
        raise NotImplementedError

    @property
    @abstractmethod
    def attribute_values(self) -> Iterable[GtfAttributeValueType]:
        """Other attributes presented in Key-Value pair"""
        raise NotImplementedError

    @property
    @abstractmethod
    def naive_length(self) -> int:
        raise NotImplementedError


@total_ordering
class Feature(FeatureInterface):
    """
    A general GTF/GFF/BED Record.
    """

    __slots__ = (
        '_seqname',
        '_source',
        '_feature',
        '_start',
        '_end',
        '_score',
        '_strand',
        '_frame',
        '_attribute',
        '_parsed_feature'
    )

    _seqname: str
    _source: Optional[str]
    _feature: Optional[str]
    _parsed_feature: Optional[FeatureType]
    _start: int
    _end: int
    _score: Optional[Union[int, float]]
    _strand: Optional[bool]
    _frame: Optional[int]
    _attribute: GTFAttributeType

    @property
    def seqname(self) -> str:
        return self._seqname

    @property
    def source(self) -> Optional[str]:
        return self._source

    @property
    def feature(self) -> Optional[str]:
        return self._feature

    @property
    def parsed_feature(self) -> FeatureType:
        if self._parsed_feature is None:
            if self._feature is None:
                self._parsed_feature = FeatureType.NOT_PRESENT
            else:
                self._parsed_feature = _raw_feature_type_translator.get(
                    self._feature.lower(),
                    FeatureType.UNKNOWN
                )
        return self._parsed_feature

    @property
    def start(self) -> int:
        return self._start

    @property
    def start0b(self) -> int:
        return self._start - 1

    @property
    def end(self) -> int:
        return self._end

    @property
    def end0b(self) -> int:
        return self._end

    @property
    def score(self) -> Optional[Union[int, float]]:
        return self._score

    @property
    def strand(self) -> Optional[bool]:
        return self._strand

    @property
    def frame(self) -> Optional[int]:
        return self._frame

    @property
    def attribute_keys(self) -> Iterable[str]:
        return self._attribute.keys()

    @property
    def attribute_values(self) -> Iterable[GtfAttributeValueType]:
        return self._attribute.values()

    @property
    def naive_length(self) -> int:
        return self.end0b - self.start0b

    def attribute_get(self, name: str, default: Optional[GtfAttributeValueType] = None) -> GtfAttributeValueType:
        """Other attributes presented in Key-Value pair"""
        return self._attribute.get(name, default)

    def overlaps(self, other: Feature) -> bool:
        if self.seqname != other.seqname:
            return False
        return (
                self.start < other.start < self.end or
                self.start < other.end < self.end or
                (
                        other.start < self.start and
                        self.end < other.end
                )
        )

    def update(
            self,
            *,
            seqname: Union[str, _NotSet] = _notset,
            source: Union[Optional[str], _NotSet] = _notset,
            feature: Union[Optional[str], _NotSet] = _notset,
            start: Union[int, _NotSet] = _notset,
            end: Union[int, _NotSet] = _notset,
            score: Union[Optional[Union[int, float]], _NotSet] = _notset,
            strand: Union[Optional[bool], _NotSet] = _notset,
            frame: Union[Optional[int], _NotSet] = _notset,
            attribute: Union[GTFAttributeType, _NotSet] = _notset
    ) -> Feature:
        return Feature(
            seqname=self._seqname if seqname is _notset else seqname,
            source=self._source if source is _notset else source,
            feature=self._feature if feature is _notset else feature,
            start=self._start if start is _notset else start,
            end=self._end if end is _notset else end,
            score=self._score if score is _notset else score,
            strand=self._strand if strand is _notset else strand,
            frame=self._frame if frame is _notset else frame,
            attribute=self._attribute if attribute is _notset else attribute,
        )

    def update_attribute(
            self,
            **attribute
    ):
        new_attribute = dict(self._attribute)
        new_attribute.update(attribute)
        return Feature(
            seqname=self._seqname,
            source=self._source,
            feature=self._feature,
            start=self._start,
            end=self._end,
            score=self._score,
            strand=self._strand,
            frame=self._frame,
            attribute=new_attribute,
        )

    def reset_attribute(self, **attribute) -> Feature:
        return Feature(
            seqname=self._seqname,
            source=self._source,
            feature=self._feature,
            start=self._start,
            end=self._end,
            score=self._score,
            strand=self._strand,
            frame=self._frame,
            attribute=attribute,
        )

    def keep_only_selected_attribute(self, *attribute_names) -> Feature:

        return Feature(
            seqname=self._seqname,
            source=self._source,
            feature=self._feature,
            start=self._start,
            end=self._end,
            score=self._score,
            strand=self._strand,
            frame=self._frame,
            attribute={k: self._attribute[k] for k in self._attribute.keys() if k in attribute_names},
        )

    def __init__(
            self,
            seqname: str,
            source: Optional[str],
            feature: Optional[str],
            start: int,
            end: int,
            score: Optional[Union[int, float]],
            strand: Optional[bool],
            frame: Optional[int],
            attribute: Optional[GTFAttributeType] = None
    ):
        """
        The filenames are named after Ensembl specifications.

        .. warning::
            Ensembl uses different way to represent 5'UTR.
        """
        if start < 1:
            raise RegionError(f"Start ({start}) cannot less than 1")
        if end < 1:
            raise RegionError(f"End ({end}) cannot less than 1")
        if end < start:
            raise RegionError(f"End ({end}) cannot less than Start ({start})")
        self._seqname = seqname
        self._source = source
        self._feature = feature
        self._parsed_feature = None
        self._start = start
        self._end = end
        self._score = score
        self._strand = strand
        self._frame = frame
        if attribute is None:
            attribute = {}
        self._attribute = dict(attribute)

    def __eq__(self, other: Feature):
        return self.seqname == other.seqname and \
               self.source == other.source and \
               self.feature == other.feature and \
               self.start == other.start and \
               self.end == other.end and \
               self.score == other.score and \
               self.strand == other.strand and \
               self.frame == other.frame and \
               list(self.attribute_keys) == list(other.attribute_keys) and \
               list(self.attribute_values) == list(other.attribute_values)

    def regional_equiv(self, other: Feature):
        return self._start == other.start and \
               self._end == other.end and \
               self._seqname == other.seqname and \
               self._strand == other.strand

    def __gt__(self, other: Feature):
        if not self.seqname == other.seqname:
            return self.seqname > other.seqname
        if not self.start == other.start:
            return self.start > other.start
        if not self.parsed_feature == other.parsed_feature:
            return self.parsed_feature > other.parsed_feature
        if not self.end == other.end:
            return self.end > other.end

    def to_dict(self) -> Dict[str, str]:
        return {
            "seqname": self.seqname,
            "source": self.source,
            "feature": self.feature,
            "start": str(self.start),
            "end": str(self.end),
            "score": str(self.score),
            "strand": str(self.strand),
            "frame": str(self.frame),
            "attribute": self._attribute
        }

    def __repr__(self):
        return str(self.to_dict())

