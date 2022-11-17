from __future__ import annotations

import bisect
from typing import List, Iterable, Optional, Tuple

from labw_utils.bioutils.datastructure.gv import GVPError, ContainerInterface
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.datastructure.gv.gene import Gene
from labw_utils.bioutils.datastructure.gv.gene_container import GeneContainerInterface
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container import TranscriptContainerInterface
from labw_utils.bioutils.record.feature import Feature, FeatureType


class DuplicatedGeneIDError(GVPError):
    pass


class GeneTree(
    GeneContainerInterface,
    TranscriptContainerInterface,
    ContainerInterface
):
    __slots__ = (
        "_genes",
        "_gene_ids"
    )
    _genes: List[Gene]
    _gene_ids: List[str]
    _transcript_ids_to_gene_ids_index:Tuple[List[str], List[str]]

    @property
    def number_of_genes(self) -> int:
        return len(self._genes)

    @property
    def gene_values(self) -> Iterable[Gene]:
        return iter(self._genes)

    @property
    def gene_ids(self) -> Iterable[str]:
        return iter(self._gene_ids)

    @property
    def number_of_transcripts(self) -> int:
        return sum(transcript_container.number_of_transcripts for transcript_container in self._genes)

    @property
    def transcript_values(self) -> Iterable[Transcript]:
        for gene in self._genes:
            yield from gene.transcript_values

    @property
    def transcript_ids(self) -> Iterable[str]:
        return iter(self._transcript_ids_to_gene_ids_index[1])

    def __init__(
            self,
            *,
            keep_sorted: bool = False,
            genes: Optional[List[Gene]] = None,
            gene_ids: Optional[List[str]] = None,
            transcript_ids_to_gene_ids_index:Optional[Tuple[List[str], List[str]]] = None
    ):
        self._is_sorted = keep_sorted
        if genes is None:
            self._genes = []
        else:
            self._genes = list(genes)
        if gene_ids is None:
            self._gene_ids = list(gene.gene_id for gene in self._genes)
        else:
            self._gene_ids = list(gene_ids)
        if transcript_ids_to_gene_ids_index is None:
            if self._genes:
                raise TypeError("Users may not set transcript_ids_to_gene_ids_index")
            self._transcript_ids_to_gene_ids_index = [], []
        else:
            self._transcript_ids_to_gene_ids_index = transcript_ids_to_gene_ids_index

    def get_gene(self, gene_id: str) -> Gene:
        return self._genes[self._gene_ids.index(gene_id)]

    def add_gene(self, gene: Gene) -> GeneTree:
        new_genes = list(self._genes)
        new_gene_ids = list(self._gene_ids)
        new_transcript_ids_to_gene_ids_index = tuple(self._transcript_ids_to_gene_ids_index)
        if gene.gene_id in self._gene_ids:
            raise DuplicatedGeneIDError(
                f"Gene ID {gene.gene_id} duplicated"
            )
        if self._is_sorted:
            new_pos = bisect.bisect_left(self._genes, gene)
            new_genes.insert(new_pos, gene)
            new_gene_ids.insert(new_pos, gene.gene_id)
        else:
            new_genes.append(gene)
            new_gene_ids.append(gene.gene_id)
        for transcript_id in gene.transcript_ids:
            new_transcript_ids_to_gene_ids_index[0].append(gene.gene_id)
            new_transcript_ids_to_gene_ids_index[1].append(transcript_id)
        return GeneTree(
            keep_sorted=self._is_sorted,
            genes=new_genes,
            gene_ids=new_gene_ids,
            transcript_ids_to_gene_ids_index=new_transcript_ids_to_gene_ids_index
        )

    def del_gene(self, gene_id: str) -> GeneTree:
        new_genes = list(self._genes)
        new_gene_ids = list(self._gene_ids)
        new_transcript_ids_to_gene_ids_index = self._transcript_ids_to_gene_ids_index
        pop_index = self._gene_ids.index(gene_id)
        _ = new_genes.pop(pop_index)
        _ = new_gene_ids.pop(pop_index)
        while True:
            try:
                pop_index = new_transcript_ids_to_gene_ids_index[0].index(gene_id)
                _ = new_transcript_ids_to_gene_ids_index[0].pop(pop_index)
                _ = new_transcript_ids_to_gene_ids_index[1].pop(pop_index)
            except ValueError:
                break
        return GeneTree(
            keep_sorted=self._is_sorted,
            genes=new_genes,
            gene_ids=new_gene_ids,
            transcript_ids_to_gene_ids_index=self._transcript_ids_to_gene_ids_index
        )

    def replace_gene(self, new_gene: Gene) -> GeneTree:
        return self.del_gene(new_gene.gene_id).add_gene(new_gene)

    def get_transcript(self, transcript_id: str) -> Transcript:
        for transcript_container in self._genes:
            try:
                return transcript_container.get_transcript(transcript_id)
            except ValueError:
                pass
        raise ValueError(f"{transcript_id} not found!")

    def add_transcript(self, transcript: Transcript) -> GeneTree:
        if transcript.gene_id in self._gene_ids:
            return self.replace_gene(
                self.get_gene(transcript.gene_id).add_transcript(transcript)
            )
        else:
            return self.add_gene(
                BaseFeatureProxy.cast_to(transcript, Gene)
            ).add_transcript(transcript)

    def del_transcript(self, transcript_id: str) -> GeneTree:
        return self.replace_gene(
            self.get_gene(self.get_transcript(transcript_id).gene_id).del_transcript(transcript_id)
        )

    def replace_transcript(self, new_transcript: Transcript) -> GeneTree:
        return self.del_transcript(new_transcript.transcript_id).add_transcript(new_transcript)

    def add_exon(self, exon: Exon) -> GeneTree:
        if exon.transcript_id in self.transcript_ids:
            return self.replace_transcript(
                self.get_transcript(exon.transcript_id).add_exon(exon)
            )
        else:
            return self.add_transcript(
                Transcript.infer_from_exon(exon)
            )

    def del_exon(self, transcript_id: str, exon_number: int) -> GeneTree:
        return self.replace_transcript(self.get_transcript(transcript_id).del_exon(exon_number))

    def replace_exon(self, transcript_id: str, exon: Exon) -> GeneTree:
        return self.replace_transcript(self.get_transcript(transcript_id).replace_exon(exon))

    def _add(self, feature: Feature) -> GeneTree:
        if feature.parsed_feature == FeatureType.Exon:
            return self.add_exon(Exon(feature))
        elif feature.parsed_feature == FeatureType.Transcript:
            return self.add_transcript(Transcript(feature))
        elif feature.parsed_feature == FeatureType.Gene:
            return self.add_gene(Gene(feature))
        else:
            return self

    @classmethod
    def from_feature_iterator(cls, feature_iterator: Iterable[Feature]):
        new_instance = cls()
        for feature in feature_iterator:
            new_instance = new_instance._add(feature)
        return new_instance
