from __future__ import annotations

import bisect
from typing import List, Iterable

from labw_utils.bioutils.datastructure.gv import GVPError
from labw_utils.bioutils.datastructure.gv.gene import Gene
from labw_utils.bioutils.datastructure.gv.gene_container import GeneContainerInterface
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container import TranscriptContainerInterface


class DuplicatedGeneIDError(GVPError):
    pass


class GeneTree(GeneContainerInterface, TranscriptContainerInterface):
    __slots__ = (
        "_genes",
        "_gene_ids"
    )
    _genes: List[Gene]
    _gene_ids: List[str]

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
        for transcript_container in self._genes:
            yield from transcript_container.transcript_values

    @property
    def transcript_ids(self) -> Iterable[str]:
        for transcript_container in self._genes:
            yield from transcript_container.transcript_ids

    def __init__(
            self,
            genes: Iterable[Gene]
    ):
        self._genes = []
        self._gene_ids = []
        self._genes = list(genes)
        self._gene_ids = list(gene.gene_id for gene in genes)

    def get_gene(self, gene_id: str) -> Gene:
        return self._genes[self._gene_ids.index(gene_id)]

    def add_gene(self, gene: Gene) -> GeneTree:
        new_genes = list(self._genes)
        if gene.gene_id in self._genes:
            raise DuplicatedGeneIDError(
                f"Gene ID {gene.gene_id} duplicated"
            )
        new_pos = bisect.bisect_left(self._genes, gene)
        new_genes.insert(new_pos, gene)
        return GeneTree(new_genes)

    def del_gene(self, gene_id: str) -> GeneTree:
        new_genes = list(self._genes)
        _ = new_genes.pop(self._gene_ids.index(gene_id))
        return GeneTree(new_genes)

    def replace_gene(self, new_gene: Gene) -> GeneTree:
        return self.del_gene(new_gene.gene_id).add_gene(new_gene)

    def get_transcript(self, transcript_id: str) -> Transcript:
        for transcript_container in self._genes:
            try:
                return transcript_container.get_transcript(transcript_id)
            except KeyError:
                pass
        raise KeyError(f"{transcript_id} not found!")

    def add_transcript(self, transcript: Transcript) -> GeneTree:
        if transcript.gene_id in self._gene_ids:
            return self.replace_gene(
                self.get_gene(transcript.gene_id).add_transcript(transcript)
            )
        else:
            pass

    def del_transcript(self, transcript_id: str) -> GeneTree:
        pass

    def replace_transcript(self, new_transcript: Transcript) -> GeneTree:
        return self.del_transcript(new_transcript.gene_id).add_transcript(new_transcript)
