"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

from labw_utils.bioutils.datastructure.region_indexer import NumpyIntervalEngine
from labw_utils.bioutils.record.feature import Feature
from labw_utils.typing_importer import Iterable, List, Tuple, Optional


class QuantificationOptimizedGeneTree:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    _feature_ids: List[str]
    _feature_boundary: NumpyIntervalEngine
    _chromosome_names: List[str]

    def iter_chromosome_names(self) -> Iterable[str]:
        return iter(self._chromosome_names)

    def __init__(self, feature_ids: List[str], feature_boundary: NumpyIntervalEngine, chromosome_names: List[str]):
        self._feature_ids = list(feature_ids)
        self._feature_boundary = feature_boundary
        self._chromosome_names = list(chromosome_names)

    @classmethod
    def from_feature_iterator(
        cls,
        feature_iterator: Iterable[Feature],
        feature_attribute_name: str = "transcript_id",
        feature_type: str = "exon",
    ) -> QuantificationOptimizedGeneTree:
        staged_features = []
        chromosome_names = set()
        for feature in feature_iterator:
            if feature.feature == feature_type and feature.attribute_get(feature_attribute_name) is not None:
                staged_features.append(feature)
                chromosome_names.add(feature.seqname)

        nie = NumpyIntervalEngine.from_interval_iterator(
            ((_feature.seqname, _feature.strand), _feature.start0b, _feature.end0b) for _feature in staged_features
        )
        feature_ids = list(_feature.attribute_get(feature_attribute_name) for _feature in staged_features)
        return cls(feature_ids, nie, list(chromosome_names))

    def overlap(self, query_interval: Tuple[Tuple[str, Optional[bool]], int, int]) -> List[str]:
        return [self._feature_ids[i] for i in self._feature_boundary.overlap(query_interval)]
