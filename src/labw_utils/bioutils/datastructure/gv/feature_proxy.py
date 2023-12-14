"""
TODO: docs

.. versionadded:: 1.0.2
"""

from __future__ import annotations

__all__ = (
    "BaseFeatureProxy",
    "update_gene_id",
    "update_transcript_id",
)

from labw_utils.bioutils.datastructure.gv import (
    CanCheckInterface,
    generate_unknown_gene_id,
    generate_unknown_transcript_id,
    SequenceFuncType,
    LegalizeRegionFuncType,
    dumb_legalize_region_func,
)
from labw_utils.bioutils.record.feature import (
    FeatureType,
    GtfAttributeValueType,
    FeatureInterface,
    BiologicalIntervalInterface,
    NotSet,
    GtfAttributeType,
    notset,
)
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Optional, TypeVar, Union, Callable, Type, Tuple, Literal, Dict, Any

from labw_utils.typing_importer import SequenceProxy

lh = get_logger(__name__)

_T = TypeVar("_T")
_OutType = TypeVar("_OutType")


class BaseFeatureProxy(FeatureInterface, CanCheckInterface):
    """
    Base class of Feature Proxy.

    .. versionadded:: 1.0.2
    """

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
    ) -> BaseFeatureProxy:
        """
        Update the data inside.

        .. warning:: This function is in-place.
        """
        self._data = self._data.update(
            seqname=seqname,
            source=source,
            feature=feature,
            start=start,
            end=end,
            score=score,
            strand=strand,
            frame=frame,
            attribute=attribute,
        )
        return self

    def update_attribute(self, **attribute) -> BaseFeatureProxy:
        """
        Update the data inside.

        .. warning:: This function is in-place.
        """
        self._data = self._data.update_attribute(**attribute)
        return self

    def reset_attribute(self, **attribute) -> BaseFeatureProxy:
        """
        Update the data inside.

        .. warning:: This function is in-place.
        """
        self._data = self._data.reset_attribute(**attribute)
        return self

    def keep_only_selected_attribute(self, *attribute_names) -> BaseFeatureProxy:
        """
        Update the data inside.

        .. warning:: This function is in-place.
        """
        self._data = self._data.keep_only_selected_attribute(**attribute_names)
        return self

    __slots__ = ("_data",)
    _data: FeatureInterface

    def __init__(self, *, data: FeatureInterface, is_checked: bool, **kwargs):
        """
        Modification of `data` will cause errors!
        """
        _ = kwargs
        self._is_checked = is_checked
        self._data = data

    def __repr__(self):
        return "BaseGeneViewFeature"

    def __str__(self):
        return repr(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BaseFeatureProxy):
            raise TypeError
        return self._data == other._data

    def __gt__(self, other: FeatureInterface) -> bool:
        return self._data > other

    def overlaps(self, other: BiologicalIntervalInterface, is_stranded: bool = True) -> bool:
        return self._data.overlaps(other, is_stranded=is_stranded)

    @property
    def naive_length(self) -> int:
        return self._data.naive_length

    @property
    def start(self) -> int:
        return self._data.start

    @property
    def start0b(self) -> int:
        return self._data.start0b

    @property
    def end(self) -> int:
        return self._data.end

    @property
    def end0b(self) -> int:
        return self._data.end0b

    @property
    def feature(self) -> Optional[str]:
        return self._data.feature

    @property
    def parsed_feature(self) -> FeatureType:
        return self._data.parsed_feature

    @property
    def seqname(self) -> str:
        return self._data.seqname

    @property
    def strand(self) -> Optional[bool]:
        return self._data.strand

    @property
    def source(self) -> Optional[str]:
        return self._data.source

    @property
    def score(self) -> Optional[Union[int, float]]:
        return self._data.score

    @property
    def frame(self) -> Optional[int]:
        return self._data.frame

    @property
    def attribute_values(self) -> SequenceProxy[GtfAttributeValueType]:
        return SequenceProxy(self._data.attribute_values)

    @property
    def attribute_keys(self) -> SequenceProxy[str]:
        return SequenceProxy(self._data.attribute_keys)

    def attribute_get(self, name: str, default: Optional[GtfAttributeValueType] = None) -> GtfAttributeValueType:
        return self._data.attribute_get(name, default)

    def attribute_get_coerce(
        self,
        name: str,
        out_type: Type[_OutType],
        coerce_func: Optional[Callable[[Optional[GtfAttributeValueType]], _OutType]] = None,
        default: Optional[GtfAttributeValueType] = None,
    ) -> _OutType:
        return self._data.attribute_get_coerce(name, out_type, coerce_func, default)

    def get_data(self) -> FeatureInterface:
        return self._data

    def regional_equiv(self, other: BiologicalIntervalInterface, is_stranded: bool = True):
        return self._data.regional_equiv(other, is_stranded)

    def flanking_sequence(
        self,
        direction: Literal["left", "right"],
        length: int,
        sequence_func: SequenceFuncType,
        legalize_region_func: LegalizeRegionFuncType = dumb_legalize_region_func,
    ) -> str:
        if direction == "left":
            return sequence_func(*legalize_region_func(self.seqname, self.start0b - length, self.start0b))
        else:
            return sequence_func(*legalize_region_func(self.seqname, self.end0b, self.end0b + length))

    def gc(self):
        """
        Perform garbage collection on cached cDNA, etc.
        """
        pass


def update_gene_id(data: FeatureInterface) -> Tuple[str, FeatureInterface]:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    gene_id = data.attribute_get("gene_id")
    if gene_id is None:
        gene_id = generate_unknown_gene_id()
        data = data.update_attribute(gene_id=gene_id)
    elif not isinstance(gene_id, str):
        gene_id = data.attribute_get_coerce("gene_id", str, str)
        data = data.update_attribute(gene_id=gene_id)
    return gene_id, data


def update_transcript_id(data: FeatureInterface) -> Tuple[str, FeatureInterface]:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    transcript_id = data.attribute_get("transcript_id")
    if transcript_id is None:
        transcript_id = generate_unknown_transcript_id()
        data = data.update_attribute(transcript_id=transcript_id)
    elif not isinstance(transcript_id, str):
        transcript_id = data.attribute_get_coerce("transcript_id", str, str)
        data = data.update_attribute(transcript_id=transcript_id)
    return transcript_id, data
