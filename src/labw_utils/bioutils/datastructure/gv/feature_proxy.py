"""
gene_view_proy -- GTF/GFF3/BED Record Proxy for Features in GeneView without Data Loss
"""

from __future__ import annotations

__all__ = [
    'TranscriptInAGeneOnDifferentStrandError',
    'DuplicatedTranscriptIDError',
    'DuplicatedTranscriptError',
    'TranscriptInAGeneOnDifferentChromosomeError'
]

import copy
from typing import Optional, Iterable, Type

from labw_utils.bioutils.datastructure._gv_errors import _all as _gve_all
from labw_utils.bioutils.datastructure.gv import _T, GVPError
from labw_utils.bioutils.record.feature import Feature, FeatureType, GtfAttributeValueType
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.devutils.decorators import copy_doc

__all__.extend(_gve_all)

lh = get_logger(__name__)


class TranscriptInAGeneOnDifferentChromosomeError(GVPError):
    pass


class DuplicatedTranscriptError(GVPError):
    pass


class DuplicatedTranscriptIDError(GVPError):
    pass


class TranscriptInAGeneOnDifferentStrandError(GVPError):
    pass


class BaseFeatureProxy:
    """
    Base class of Feature Proxy.
    """
    __slots__ = (
        "_data",
    )
    _data: Feature
    _preserved_attrs: Iterable[str]

    def cast_to(
            self,
            class_type: Type[_T],
            **kwargs
    ) -> _T:
        """
        Cast the type of current feature proxy to another using existing data.

        :param class_type: Destination type.
        """
        return class_type(self._data, **kwargs)

    def __init__(self, data: Feature, **kwargs):
        if not hasattr(self, "_preserved_attrs"):
            # TODO
            raise TypeError("_preserved_attrs not defined; ask your maintainer for this problem")
        self._data = copy.deepcopy(data)

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
        data_attributes_update_kwargs = {k: v for k, v in zip(
            self._data.attribute_keys, self._data.attribute_values
        )}
        for preserved_attribute_name in self._preserved_attrs:
            if preserved_attribute_name in kwargs:
                data_attributes_update_kwargs[preserved_attribute_name] = kwargs[preserved_attribute_name]
        new_data = self._data.update(attribute=data_attributes_update_kwargs)
        return self.__class__(new_data, shortcut=True)

    def update_data(self, **kwargs) -> BaseFeatureProxy:
        if "attribute" in kwargs and any(
                map(
                    lambda preserved_attribute_name: preserved_attribute_name in kwargs["attribute"],
                    self._preserved_attrs
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

#
# class Gene(BaseFeatureProxy):
#     _transcript_ids: List[str]
#     _transcripts: List[Transcript]
#
#     @property
#     def transcribed_length(self):
#         raise TypeError("Do not know how to get transcribed length for a gene")
#
#     def check_transcript_duplication(self) -> Optional[Tuple[str, str]]:
#         for i in range(self.number_of_transcripts):
#             for j in range(i, self.number_of_transcripts):
#                 if self._transcripts[i] == self._transcripts[j]:
#                     return self._transcripts[i].transcript_id, self._transcripts[j].transcript_id
#         return None
#
#     def check_whether_one_transcript_duplicates_with_others(self, transcript_id: str) -> Optional[str]:
#         transcript = self.get_transcript(transcript_id)
#         for other_transcript in self.iter_transcripts():
#             if other_transcript == transcript and \
#                     other_transcript.transcript_id != transcript.transcript_id:
#                 return other_transcript.transcript_id
#         return None
#
#     @property
#     def gene_id(self) -> str:
#         return self._data.attribute["gene_id"]
#
#     @gene_id.setter
#     def gene_id(self, value: str):
#         self._data.attribute["gene_id"] = value
#
#     @property
#     def number_of_transcripts(self) -> int:
#         return len(self._transcripts)
#
#     def get_transcript(self, transcript_id: str) -> Transcript:
#         return self._transcripts[self._transcript_ids.index(transcript_id)]
#
#     def iter_transcripts(self) -> Iterable[Transcript]:
#         return self._transcripts
#
#     def iter_transcript_ids(self) -> Iterable[str]:
#         return self._transcript_ids
#
#     def _setup(self):
#         self._transcripts = []
#         self._transcript_ids = []
#
#     def _setup_gtf(self) -> None:
#         if "gene_id" not in self._data.attribute:
#             self._data.attribute["gene_id"] = generate_unknown_gene_id()
#
#     def _setup_gff3(self) -> None:
#         raise NotImplementedError
#
#     def __repr__(self):
#         return f"Gene {self.gene_id}"
