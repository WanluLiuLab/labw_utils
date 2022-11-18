"""

>>> import numpy as np
>>> splice_sites = [[1, 2], [4, 5], [6, 7]]
>>> splice_sites_np = np.array(splice_sites, dtype=int)
>>> splice_sites_np
array([[1, 2],
       [4, 5],
       [6, 7]])
>>> splice_sites_np.shape
(3, 2)
"""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, List, Tuple, Dict

import numpy as np
import numpy.typing as npt

from labw_utils.bioutils.record.feature import Feature, FeatureType


class QuantificationOptimizedGene:
    _transcript_ids: List[str]
    _exon_boundaries: npt.NDArray[int]
    _exon_boundaries_graph: npt.NDArray[bool]
    _splice_sites: npt.NDArray[int]
    _splice_sites_graph: npt.NDArray[bool]
    _self_boundaries: Tuple[int, int]

    @classmethod
    def from_exons(cls, exons: List[Feature]) -> QuantificationOptimizedGene:
        transcript_id_to_exon_index = defaultdict(lambda: [])
        for exon in exons:
            transcript_id_to_exon_index[exon.attribute_get("transcript_id")].append(exon)

        transcript_id_to_splice_site_index: Dict[str, List[Tuple[int, int]]] = {}
        transcript_id_to_exon_boundary_index: Dict[str, List[Tuple[int, int]]] = {}
        for transcript_id, selected_exons in transcript_id_to_exon_index.items():
            selected_exons = sorted(selected_exons)
            transcript_id_to_splice_site_index[transcript_id] = list(
                (selected_exons[i].end, selected_exons[i + 1].start)
                for i in range(len(selected_exons) - 1)
            )
            transcript_id_to_exon_boundary_index[transcript_id] = list(
                (exon.start, exon.end) for exon in selected_exons
            )
        splice_sites: List[Tuple[int, int]] = list(set(*transcript_id_to_splice_site_index.values()))
        exon_boundaries: List[Tuple[int, int]] = list(set(*transcript_id_to_exon_boundary_index.values()))
        transcript_ids: List[str] = list(set(transcript_id_to_exon_index.keys()))
        splice_sites_graph = np.ndarray((len(splice_sites), len(transcript_ids)), dtype=bool)
        exon_boundaries_graph = np.ndarray((len(exon_boundaries), len(transcript_ids)), dtype=bool)


class QuantificationOptimizedGeneTree:
    _gene_ids: List[str]
    _transcript_boundary: npt.NDArray

    @classmethod
    def from_feature_iterator(
            cls,
            feature_iterator: Iterable[Feature]
    ) -> QuantificationOptimizedGeneTree:
        gene_id_to_exon_index = defaultdict(lambda: [])
        """gene_id -> List[Exon]"""
        for feature in feature_iterator:
            if feature.parsed_feature != FeatureType.Exon:
                continue
            feature = feature.keep_only_selected_attribute("gene_id", "transcript_id")
            if feature.attribute_get("gene_id") is None or feature.attribute_get("transcript_id") is None:
                raise ValueError("Feature needs gene_id and transcript_id")
            gene_id_to_exon_index[feature.attribute_get("gene_id")].append(feature)

        for exon in gene_id_to_exon_index:
            exon.attribute_get("gene_id")

        ...
