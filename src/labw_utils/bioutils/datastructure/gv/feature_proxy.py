"""
gene_view_proy -- GTF/GFF3/BED Record Proxy for Features in GeneView without Data Loss
"""

from __future__ import annotations

__all__ = [
    'BaseFeatureProxy'
]

from typing import Optional, Iterable, TypeVar, Union

from labw_utils.bioutils.datastructure.gv import CanCheckInterface
from labw_utils.bioutils.record.feature import Feature, FeatureType, GtfAttributeValueType, FeatureInterface
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)

_T = TypeVar("_T")


class BaseFeatureProxy(FeatureInterface, CanCheckInterface):
    """
    Base class of Feature Proxy.
    """

    __slots__ = [
        "_data",
    ]
    _data: Feature

    def __init__(
            self,
            *,
            data: Feature,
            is_checked: bool,
            **kwargs
    ):
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

    def __eq__(self, other: BaseFeatureProxy) -> bool:
        return self._data == other._data

    def __ne__(self, other: BaseFeatureProxy) -> bool:
        return self._data != other._data

    def __gt__(self, other: BaseFeatureProxy) -> bool:
        return self._data > other._data

    def __ge__(self, other: BaseFeatureProxy) -> bool:
        return self._data >= other._data

    def __lt__(self, other: BaseFeatureProxy) -> bool:
        return self._data < other._data

    def __le__(self, other: BaseFeatureProxy) -> bool:
        return self._data <= other._data

    def overlaps(self, other: BaseFeatureProxy) -> bool:
        return self._data.overlaps(other._data)

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
    def feature(self) -> str:
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
    def attribute_values(self) -> Iterable[GtfAttributeValueType]:
        return self._data.attribute_values

    @property
    def attribute_keys(self) -> Iterable[str]:
        return self._data.attribute_keys

    def attribute_get(self, name: str, default: Optional[GtfAttributeValueType] = None) -> GtfAttributeValueType:
        return self._data.attribute_get(name, default)

    def get_data(self) -> Feature:
        return self._data
