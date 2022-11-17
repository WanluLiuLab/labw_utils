"""
gene_view_proy -- GTF/GFF3/BED Record Proxy for Features in GeneView without Data Loss
"""

from __future__ import annotations

__all__ = [
    'BaseFeatureProxy'
]

from typing import Optional, Iterable, Type, TypeVar

from labw_utils.bioutils.record.feature import Feature, FeatureType, GtfAttributeValueType
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.devutils.decorators import copy_doc

lh = get_logger(__name__)

_T = TypeVar("_T")


class BaseFeatureProxy:
    """
    Base class of Feature Proxy.
    """

    __slots__ = (
        "_data",
    )
    _data: Feature
    preserved_attributes: Iterable[str]

    def cast_to(
            self,
            class_type: Type[_T],
            **kwargs
    ) -> _T:
        """
        Cast the type of current feature proxy to another using existing data.

        :param class_type: Destination type.

        Example:

            BaseFeatureProxy.cast_to(exon, class_type=Transcript, is_inferred=True)

        would cast an exon to a transcript.
        """
        return class_type(self._data, **kwargs)

    def __init__(self, data: Feature, **kwargs):
        """
        Modification of `data` will cause errors!
        """
        if not hasattr(self, "preserved_attributes"):
            raise TypeError("preserved_attributes not defined; ask your maintainer for this problem")
        self._data = data

    def __repr__(self):
        return "BaseGeneViewFeature"

    def __str__(self):
        return repr(self)

    @copy_doc(Feature.__eq__)
    def __eq__(self, other: BaseFeatureProxy) -> bool:
        return self._data == other._data

    @copy_doc(Feature.__ne__)
    def __ne__(self, other: BaseFeatureProxy) -> bool:
        return self._data != other._data

    @copy_doc(Feature.__gt__)
    def __gt__(self, other: BaseFeatureProxy) -> bool:
        return self._data > other._data

    @copy_doc(Feature.__ge__)
    def __ge__(self, other: BaseFeatureProxy) -> bool:
        return self._data >= other._data

    @copy_doc(Feature.__lt__)
    def __lt__(self, other: BaseFeatureProxy) -> bool:
        return self._data < other._data

    @copy_doc(Feature.__le__)
    def __le__(self, other: BaseFeatureProxy) -> bool:
        return self._data <= other._data

    def update(self, **kwargs) -> BaseFeatureProxy:
        """
        Update data attributes that ARE related to ``preserved_attributes``.

        :param kwargs: Arguments for :py:func:`Feature.update`
        """
        data_attributes_update_kwargs = {k: v for k, v in zip(
            self._data.attribute_keys, self._data.attribute_values
        )}
        for preserved_attribute_name in self.preserved_attributes:
            if preserved_attribute_name in kwargs:
                data_attributes_update_kwargs[preserved_attribute_name] = kwargs[preserved_attribute_name]
        new_data = self._data.update(attribute=data_attributes_update_kwargs)
        return self.__class__(new_data, shortcut=True)

    def update_data(self, **kwargs) -> BaseFeatureProxy:
        """
        Update data attributes that are NOT related to ``preserved_attributes``.

        :param kwargs: Arguments for :py:func:`Feature.update`
        """
        if "attribute" in kwargs and any(
                map(
                    lambda preserved_attribute_name: preserved_attribute_name in kwargs["attribute"],
                    self.preserved_attributes
                )
        ):
            raise ValueError("Preserved attribute should not be updated here")
        new_data = self._data.update(**kwargs)
        return self.__class__(new_data, shortcut=True)

    @copy_doc(Feature.overlaps)
    def overlaps(self, other: BaseFeatureProxy) -> bool:
        return self._data.overlaps(other._data)

    @copy_doc(Feature.naive_length)
    @property
    def naive_length(self) -> int:
        return self._data.naive_length

    @copy_doc(Feature.start)
    @property
    def start(self) -> int:
        return self._data.start

    @copy_doc(Feature.start0b)
    @property
    def start0b(self) -> int:
        return self._data.start0b

    @copy_doc(Feature.end)
    @property
    def end(self) -> int:
        return self._data.end

    @copy_doc(Feature.end0b)
    @property
    def end0b(self) -> int:
        return self._data.end0b

    @copy_doc(Feature.feature)
    @property
    def feature(self) -> str:
        return self._data.feature

    @copy_doc(Feature.parsed_feature)
    @property
    def parsed_feature(self) -> FeatureType:
        return self._data.parsed_feature

    @copy_doc(Feature.seqname)
    @property
    def seqname(self) -> str:
        return self._data.seqname

    @copy_doc(Feature.strand)
    @property
    def strand(self) -> Optional[bool]:
        return self._data.strand

    @copy_doc(Feature.frame)
    @property
    def frame(self) -> Optional[int]:
        return self._data.frame

    @copy_doc(Feature.attribute_values)
    @property
    def attribute_values(self) -> Iterable[GtfAttributeValueType]:
        return self._data.attribute_values

    @copy_doc(Feature.attribute_keys)
    @property
    def attribute_keys(self) -> Iterable[str]:
        return self._data.attribute_keys

    @copy_doc(Feature.attribute_get)
    def attribute_get(self, name: str, default: Optional[GtfAttributeValueType] = None) -> GtfAttributeValueType:
        return self._data.attribute_get(name, default)
