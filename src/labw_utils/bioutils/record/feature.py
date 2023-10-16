"""
``labw_utils.bioutils.record.feature`` -- General-Purposed GTF/GFF3/BED Record that Represents a Genomic Feature

This module includes GTF/GFF3/BED record datastructure and their one-line parsers.

.. versionadded:: 1.0.2
"""

from __future__ import annotations

import enum
from abc import abstractmethod, ABC
from functools import total_ordering

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.devutils.decorators import create_class_init_doc_from_property
from labw_utils.typing_importer import Union, Optional, Mapping, List, Type, TypeVar, Callable, Sequence

from labw_utils.typing_importer import SequenceProxy


_lh = get_logger(__name__)

GtfAttributeValueType = Union[str, int, float, bool, None, List[str], List[int], List[float], List[bool], List[None]]
"""
Type of GTF/GFF attributes

.. versionadded:: 1.0.2
"""

GtfAttributeType = Mapping[str, GtfAttributeValueType]
"""
Type of GTF/GFF fields

.. versionadded:: 1.0.2
"""

VALID_GTF_QUOTE_OPTIONS = ("none", "blank", "string", "all")
"""
Valid GTF Quoting Options. They are:

* ``none`` Never quote, even if blanks found inside.
* ``blank``: Quote if blanks (\\r, \\n, \\t, space, etc.) found inside.
* ``string``: Quote if the field have type string. Will not quote for numeric types.
* ``all``: Quote all fields.

.. versionadded:: 1.0.2
"""

DEFAULT_GTF_QUOTE_OPTIONS = "all"
"""
.. versionadded:: 1.0.2
"""

_OutType = TypeVar("_OutType")


class NotSet:
    """
    Not set indicator

    .. versionadded:: 1.0.2
    """

    pass


notset = NotSet()
"""
Not set indicator

.. versionadded:: 1.0.2
"""


def feature_repr(v: GtfAttributeValueType) -> str:
    """
    Python standard :py:func:`repr` for genomic data.

    If the genomic attribute is of type ``list``, ivoke this function multiple times.

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

    .. versionadded:: 1.0.2
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


def strand_repr(strand: Optional[bool]) -> str:
    """
    Convert strand in bool to strand in str.

    .. versionadded:: 1.0.2
    """
    if strand is None:
        return "."
    if strand:
        return "+"
    return "-"


class FeatureParserError(ValueError):
    """
    General feature parsing errors.

    .. versionadded:: 1.0.2
    """

    pass


class RegionError(FeatureParserError):
    """
    Error raised if:

    1. ``start`` less than 1
    2. ``end`` less than ``start``

    .. versionadded:: 1.0.2
    """

    def __init__(self, *args):
        super(RegionError, self).__init__(*args)


class FeatureType(enum.IntEnum):
    """
    Type of genomic feature.

    .. versionadded:: 1.0.2
    """

    NOT_PRESENT = -1
    """If feature is ``.``"""

    UNKNOWN = 0
    """If not recognized."""

    EXON = 9

    FIVE_PRIME_UTR = 8
    THREE_PRIME_UTR = 7
    OTHER_UTR = 6
    """If feature is ``utr``"""

    CDS = 5
    """Coding sequence."""

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


@total_ordering
class BiologicalIntervalInterface(ABC):
    """
    Interface representing biological intervals.

    .. versionadded:: 1.0.2
    """

    @property
    @abstractmethod
    def naive_length(self) -> int:
        """Naive length, is ``self.end0b - self.start0b``"""
        raise NotImplementedError

    @property
    @abstractmethod
    def seqname(self) -> str:
        """
        Chromosome, Contig or Scaffold name.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def start(self) -> int:
        """
        Inclusive 1-based start position.
        As-is on GTF.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def start0b(self) -> int:
        """
        Inclusive 0-based start position.
        For picking from FASTA.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def end(self) -> int:
        """
        Inclusive 1-based end position.
        As-is on GTF.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def end0b(self) -> int:
        """
        Exclusive 0-based end position.
        For picking from FASTA.
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
    def regional_equiv(self, other: BiologicalIntervalInterface, is_stranded: bool = True) -> bool:
        """
        Whether two BiologicalIntervalInterface is regionally equivalent.
        Requires same ``start``, ``end`` and ``strand``.
        """
        raise NotImplementedError

    @abstractmethod
    def overlaps(self, other: BiologicalIntervalInterface, is_stranded: bool = True) -> bool:
        """
        Whether this feature overlaps with another feature.

        :param other: Another feature.
        :param is_stranded: Whether to consider strands.
        """
        raise NotImplementedError

    @abstractmethod
    def __gt__(self, other: FeatureInterface) -> bool:
        raise NotImplementedError


class BiologicalInterval(BiologicalIntervalInterface):
    """
    .. versionadded:: 1.0.3
    """

    __slots__ = (
        "_seqname",
        "_start",
        "_end",
        "_strand",
    )

    _seqname: str
    _start: int
    _end: int
    _strand: Optional[bool]

    def __init__(self, seqname: str, start: int, end: int, strand: Optional[bool]):
        if start < 1:
            raise RegionError(f"Start ({start}) cannot less than 1")
        if end < 1:
            raise RegionError(f"End ({end}) cannot less than 1")
        if end < start:
            raise RegionError(f"End ({end}) cannot less than Start ({start})")
        self._seqname = seqname
        self._start = start
        self._end = end
        self._strand = strand

    @property
    def seqname(self) -> str:
        return self._seqname

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
    def strand(self) -> Optional[bool]:
        return self._strand

    def __gt__(self, other: FeatureInterface):
        if self.seqname != other.seqname:
            return self.seqname > other.seqname
        if self.start != other.start:
            return self.start > other.start
        if self.end != other.end:
            return self.end > other.end

    def regional_equiv(self, other: BiologicalIntervalInterface, is_stranded: bool = True) -> bool:
        if is_stranded and self._strand != other.strand:
            return False
        return self._start == other.start and self._end == other.end and self._seqname == other.seqname

    def overlaps(self, other: BiologicalIntervalInterface, is_stranded: bool = True) -> bool:
        if self.seqname != other.seqname:
            return False
        if is_stranded and self.strand != other.strand:
            return False
        return (
            self.start < other.start < self.end
            or self.start < other.end < self.end
            or (other.start < self.start and self.end < other.end)
        )

    @property
    def naive_length(self) -> int:
        return self.end0b - self.start0b

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.seqname}':{self.start0b}-{self.end0b}({strand_repr(self.strand)}))"


class FeatureInterface(BiologicalIntervalInterface):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

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
        """
        Feature type that is parsed into standard format.
        """
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
        One of:

        - ``0`` (first base of the feature is the first base of a codon),
        - ``1`` (the second base is the first base of a codon),
        - ``2`` (so on).

        See Ensemble Documentations at <https://www.ensembl.org/info/website/upload/gff.html>
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def attribute_keys(self) -> Sequence[str]:
        """Other attributes presented in Key-Value pair"""
        raise NotImplementedError

    @property
    @abstractmethod
    def attribute_values(self) -> Sequence[GtfAttributeValueType]:
        """Other attributes presented in Key-Value pair"""
        raise NotImplementedError

    @abstractmethod
    def attribute_get(self, name: str, default: Optional[GtfAttributeValueType] = None) -> GtfAttributeValueType:
        """
        Get attribute with defaults and coerce.

        :param name: Attribute name.
        :param default: The default value.
        """
        raise NotImplementedError

    @abstractmethod
    def attribute_get_coerce(
        self,
        name: str,
        out_type: Type[_OutType],
        coerce_func: Optional[Callable[[Optional[GtfAttributeValueType]], _OutType]] = None,
        default: Optional[GtfAttributeValueType] = None,
    ) -> _OutType:
        """
        Get attribute with defaults and coerce.

        :param name: Attribute name.
        :param out_type: Type of output.
        :param coerce_func: Function that coerce :py:obj:`GtfAttributeValueType` to output type.
            This function should only raise :py:class:`TypeError`.
        :param default: The default value.
        :raises TypeError: On coerce failures.
        """
        raise NotImplementedError

    @abstractmethod
    def update(
        self,
        *,
        seqname: Union[str, NotSet] = notset,
        source: Union[Optional[str], NotSet] = notset,
        feature: Union[Optional[str], NotSet] = notset,
        start: Union[int, NotSet] = notset,
        end: Union[int, NotSet] = notset,
        score: Union[Optional[Union[int, float]], NotSet] = notset,
        strand: Union[Optional[bool], NotSet] = notset,
        frame: Union[Optional[int], NotSet] = notset,
        attribute: Union[GtfAttributeType, NotSet] = notset,
    ) -> FeatureInterface:
        """
        Update fields or attributes of current Feature and generate a new one.
        """
        raise NotImplementedError

    @abstractmethod
    def update_attribute(self, **attribute) -> FeatureInterface:
        """
        Perform :py:func:`dict.update` on attributes of current Feature and generate a new one.
        """
        raise NotImplementedError

    @abstractmethod
    def reset_attribute(self, **attribute) -> FeatureInterface:
        """
        Reset attributes of current Feature and generate a new one.
        """
        raise NotImplementedError

    @abstractmethod
    def keep_only_selected_attribute(self, *attribute_names) -> FeatureInterface:
        """
        Filter attributes of current Feature and generate a new one.
        """
        raise NotImplementedError


@create_class_init_doc_from_property(
    text_after="""

:raises RegionError: If the region is invalid.
"""
)
class Feature(BiologicalInterval, FeatureInterface):
    """
    A general GTF/GFF/BED Record.

    The filenames are named after Ensembl specifications.

    .. warning::
        Ensembl uses different way to represent 5'UTR.

    .. versionadded:: 1.0.2
    """

    __slots__ = (
        "_seqname",
        "_source",
        "_feature",
        "_start",
        "_end",
        "_score",
        "_strand",
        "_frame",
        "_attribute",
        "_parsed_feature",
    )

    _source: Optional[str]
    _feature: Optional[str]
    _parsed_feature: Optional[FeatureType]
    _score: Optional[Union[int, float]]
    _frame: Optional[int]
    _attribute: GtfAttributeType

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
                try:
                    self._parsed_feature = _raw_feature_type_translator[self._feature.lower()]
                except KeyError:
                    self._parsed_feature = FeatureType.UNKNOWN
        return self._parsed_feature

    @property
    def score(self) -> Optional[Union[int, float]]:
        return self._score

    @property
    def frame(self) -> Optional[int]:
        return self._frame

    @property
    def attribute_keys(self) -> Sequence[str]:
        return SequenceProxy(self._attribute.keys())

    @property
    def attribute_values(self) -> Sequence[GtfAttributeValueType]:
        return SequenceProxy(self._attribute.values())

    def attribute_get(self, name: str, default: Optional[GtfAttributeValueType] = None) -> GtfAttributeValueType:
        """Other attributes presented in Key-Value pair"""
        return self._attribute.get(name, default)

    def attribute_get_coerce(
        self,
        name: str,
        out_type: Type[_OutType],
        coerce_func: Optional[Callable[[Optional[GtfAttributeValueType]], _OutType]] = None,
        default: Optional[GtfAttributeValueType] = None,
    ) -> _OutType:
        v = self._attribute.get(name, default)
        if isinstance(v, out_type):
            return v
        if coerce_func is not None:
            return coerce_func(v)
        raise TypeError(f"Type {type(v)} cannot be coerced to type {out_type}")

    def update(
        self,
        *,
        seqname: Union[str, NotSet] = notset,
        source: Union[Optional[str], NotSet] = notset,
        feature: Union[Optional[str], NotSet] = notset,
        start: Union[int, NotSet] = notset,
        end: Union[int, NotSet] = notset,
        score: Union[Optional[Union[int, float]], NotSet] = notset,
        strand: Union[Optional[bool], NotSet] = notset,
        frame: Union[Optional[int], NotSet] = notset,
        attribute: Union[GtfAttributeType, NotSet] = notset,
    ) -> Feature:
        return Feature(
            seqname=self._seqname if isinstance(seqname, NotSet) else seqname,
            source=self._source if isinstance(source, NotSet) else source,
            feature=self._feature if isinstance(feature, NotSet) else feature,
            start=self._start if isinstance(start, NotSet) else start,
            end=self._end if isinstance(end, NotSet) else end,
            score=self._score if isinstance(score, NotSet) else score,
            strand=self._strand if isinstance(strand, NotSet) else strand,
            frame=self._frame if isinstance(frame, NotSet) else frame,
            attribute=self._attribute if isinstance(attribute, NotSet) else attribute,
        )

    def update_attribute(self, **attribute):
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
        attribute: Optional[GtfAttributeType] = None,
    ):
        super().__init__(seqname=seqname, start=start, end=end, strand=strand)

        self._source = source
        self._feature = feature
        self._parsed_feature = None
        self._score = score
        self._frame = frame
        if attribute is None:
            attribute = {}
        self._attribute = dict(attribute)

    def __eq__(self, other: object):
        if not isinstance(other, Feature):
            raise TypeError
        return (
            self.seqname == other.seqname
            and self.source == other.source
            and self.feature == other.feature
            and self.start == other.start
            and self.end == other.end
            and self.score == other.score
            and self.strand == other.strand
            and self.frame == other.frame
            and list(self.attribute_keys) == list(other.attribute_keys)
            and list(self.attribute_values) == list(other.attribute_values)
        )

    def __gt__(self, other: FeatureInterface):
        if self.seqname != other.seqname:
            return self.seqname > other.seqname
        if self.start != other.start:
            return self.start > other.start
        if self.parsed_feature != other.parsed_feature:
            return self.parsed_feature > other.parsed_feature
        if self.end != other.end:
            return self.end > other.end

    def to_dict(self) -> Mapping[str, Union[str, GtfAttributeType]]:
        """
        Convert to dict to be used for initiators.
        """
        return {
            "seqname": self.seqname,
            "source": self.source if self.source is not None else ".",
            "feature": self.feature if self.feature is not None else ".",
            "start": str(self.start),
            "end": str(self.end),
            "score": str(self.score),
            "strand": strand_repr(self.strand),
            "frame": str(self.frame) if self.frame is not None else ".",
            "attribute": self._attribute,
        }

    def __repr__(self):
        return str(self.to_dict())
