from __future__ import annotations

from typing import Iterable, List, Tuple, Optional

from labw_utils.bioutils.datastructure.region_indexer import NumpyIntervalEngine
from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.bioutils.record.feature import Feature


class QuantificationOptimizedGeneTree:
    _feature_ids: List[str]
    _feature_boundary: NumpyIntervalEngine

    def __init__(self, feature_ids: List[str], feature_boundary: NumpyIntervalEngine):
        self._feature_ids = feature_ids
        self._feature_boundary = feature_boundary

    @classmethod
    def from_feature_iterator(
            cls,
            feature_iterator: Iterable[Feature],
            feature_attribute_name: str = "transcript_id",
            feature_type: str = "exon"
    ) -> QuantificationOptimizedGeneTree:
        staged_features = []
        for feature in feature_iterator:
            if feature.feature == feature_type and feature.attribute_get(feature_attribute_name) is not None:
                staged_features.append(feature)
        nie = NumpyIntervalEngine.from_interval_iterator(
            ((_feature.seqname, _feature.strand), _feature.start0b, _feature.end0b) for _feature in staged_features
        )
        feature_ids = list(_feature.attribute_get(feature_attribute_name) for _feature in staged_features)
        return cls(feature_ids, nie)

    def overlap(self, query_interval: Tuple[Tuple[str, Optional[bool]], int, int]) -> List[str]:
        return [self._feature_ids[i] for i in self._feature_boundary.overlap(query_interval)]


if __name__ == "__main__":
    qgt = QuantificationOptimizedGeneTree.from_feature_iterator(
        GtfIterator("/home/yuzj/Documents/cpptetgs_experimental/test_data/gtf/hg38.ncbiRefSeq_sel.gtf")
    )
    print(qgt.overlap((("chr1", True), 0, 100000)))
