from __future__ import annotations

from labw_utils.bioutils.datastructure.gv import GVPError, SortedContainerInterface, CanCheckInterface
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.gene import Gene
from labw_utils.bioutils.datastructure.gv.gene_container_interface import GeneContainerInterface
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container_interface import TranscriptContainerInterface
from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.bioutils.record.feature import Feature, FeatureType
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper import pickle_helper
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.io.file_system import should_regenerate
from labw_utils.typing_importer import Iterable, Dict, Iterator

_lh = get_logger(__name__)

GVPKL_VERSION = "1.0"
"""Current version of GVPKL standard."""


class DuplicatedGeneIDError(GVPError):
    def __init__(self, gene_id: str):
        super().__init__(f"Gene ID {gene_id} duplicated")


class GeneTree(
    GeneContainerInterface,
    TranscriptContainerInterface,
    SortedContainerInterface,
    CanCheckInterface
):
    __slots__ = (
        "_gene_id_to_gene_index",
        "_transcript_ids_to_gene_ids_index"
    )
    _gene_id_to_gene_index: Dict[str, Gene]
    _transcript_ids_to_gene_ids_index: Dict[str, str]
    """transcript_id -> gene_id"""

    @property
    def number_of_genes(self) -> int:
        return len(self._gene_id_to_gene_index)

    @property
    def gene_values(self) -> Iterable[Gene]:
        return self._gene_id_to_gene_index.values()

    @property
    def gene_ids(self) -> Iterable[str]:
        return self._gene_id_to_gene_index.keys()

    @property
    def number_of_transcripts(self) -> int:
        return sum(gene.number_of_transcripts for gene in self.gene_values)

    @property
    def transcript_values(self) -> Iterable[Transcript]:
        for gene in self.gene_values:
            yield from gene.transcript_values

    @property
    def transcript_ids(self) -> Iterable[str]:
        return iter(self._transcript_ids_to_gene_ids_index.keys())

    def __init__(
            self,
            *,
            keep_sorted: bool,
            is_checked: bool,
            gene_id_to_gene_index: Dict[str, Gene],
            transcript_ids_to_gene_ids_index: Dict[str, str]
    ):
        self._is_sorted = keep_sorted
        self._is_checked = is_checked
        self._gene_id_to_gene_index = dict(gene_id_to_gene_index)
        self._transcript_ids_to_gene_ids_index = transcript_ids_to_gene_ids_index

    def get_gene(self, gene_id: str) -> Gene:
        return self._gene_id_to_gene_index[gene_id]

    def add_gene(self, gene: Gene) -> GeneTree:
        if self._is_checked and gene.gene_id in self._gene_id_to_gene_index:
            raise DuplicatedGeneIDError(gene.gene_id)
        new_gene_id_to_gene_index = dict(self._gene_id_to_gene_index)
        new_transcript_ids_to_gene_ids_index = dict(self._transcript_ids_to_gene_ids_index)
        new_gene_id_to_gene_index[gene.gene_id] = gene
        for transcript_id in gene.transcript_ids:
            new_transcript_ids_to_gene_ids_index[transcript_id] = gene.gene_id
        return GeneTree(
            keep_sorted=self._is_sorted,
            is_checked=self._is_checked,
            gene_id_to_gene_index=new_gene_id_to_gene_index,
            transcript_ids_to_gene_ids_index=new_transcript_ids_to_gene_ids_index
        )

    def del_gene(self, gene_id: str) -> GeneTree:
        new_gene_id_to_gene_index = dict(self._gene_id_to_gene_index)
        _ = new_gene_id_to_gene_index.pop(gene_id)
        new_transcript_ids_to_gene_ids_index = {
            k: v for k, v in self._transcript_ids_to_gene_ids_index.items() if v != gene_id
        }
        return GeneTree(
            keep_sorted=self._is_sorted,
            is_checked=self._is_checked,
            gene_id_to_gene_index=new_gene_id_to_gene_index,
            transcript_ids_to_gene_ids_index=new_transcript_ids_to_gene_ids_index
        )

    def replace_gene(self, new_gene: Gene) -> GeneTree:
        return self.del_gene(new_gene.gene_id).add_gene(new_gene)

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self.get_gene(
            self._transcript_ids_to_gene_ids_index[transcript_id]
        ).get_transcript(transcript_id)

    def add_transcript(self, transcript: Transcript) -> GeneTree:
        if transcript.gene_id in self._gene_id_to_gene_index:
            gene_to_be_modified = self._gene_id_to_gene_index[transcript.gene_id]
            gene_to_be_modified = gene_to_be_modified.add_transcript(transcript)
            self._gene_id_to_gene_index[gene_to_be_modified.gene_id] = gene_to_be_modified
            self._transcript_ids_to_gene_ids_index[transcript.transcript_id] = gene_to_be_modified.gene_id
            return self
        else:
            return self.add_gene(
                Gene(
                    data=transcript.get_data(),
                    is_checked=self._is_checked,
                    is_inferred=True,
                    keep_sorted=self._is_sorted,
                    shortcut=False,
                    transcripts=[],
                    transcript_ids=[]
                )
            ).add_transcript(transcript)

    def del_transcript(self, transcript_id: str) -> GeneTree:
        gene_to_be_modified = self._gene_id_to_gene_index[
            self._transcript_ids_to_gene_ids_index[transcript_id]
        ]
        gene_to_be_modified = gene_to_be_modified.del_transcript(transcript_id)
        self._gene_id_to_gene_index[gene_to_be_modified.gene_id] = gene_to_be_modified
        _ = self._transcript_ids_to_gene_ids_index.pop(transcript_id)
        return self

    def replace_transcript(self, new_transcript: Transcript) -> GeneTree:
        return self.del_transcript(new_transcript.transcript_id).add_transcript(new_transcript)

    def add_exon(self, exon: Exon) -> GeneTree:
        if exon.transcript_id in self._transcript_ids_to_gene_ids_index:
            return self.replace_transcript(
                self.get_transcript(exon.transcript_id).add_exon(exon)
            )
        else:
            return self.add_transcript(
                Transcript(
                    data=exon.get_data(),
                    exons=[],
                    is_inferred=True,
                    is_checked=self._is_checked,
                    keep_sorted=self._is_sorted,
                    shortcut=False
                )
            ).add_exon(exon)

    def del_exon(self, transcript_id: str, exon_index: int) -> GeneTree:
        return self.replace_transcript(self.get_transcript(transcript_id).del_exon(exon_index))

    def _add(self, feature: Feature) -> GeneTree:
        if feature.parsed_feature == FeatureType.EXON:
            return self.add_exon(Exon(
                data=feature,
                is_checked=self._is_checked,
                shortcut=False
            ))
        elif feature.parsed_feature == FeatureType.TRANSCRIPT:
            return self.add_transcript(Transcript(
                data=feature,
                is_checked=self._is_checked,
                keep_sorted=self._is_sorted,
                exons=[],
                is_inferred=False,
                shortcut=False
            ))
        elif feature.parsed_feature == FeatureType.GENE:
            return self.add_gene(Gene(
                data=feature,
                is_checked=self._is_checked,
                keep_sorted=self._is_sorted,
                transcripts=[],
                transcript_ids=[],
                is_inferred=False,
                shortcut=False
            ))
        else:
            return self

    @classmethod
    def from_feature_iterator_legacy(
            cls,
            feature_iterator: Iterable[Feature],
            keep_sorted: bool = False,
            is_checked: bool = False,
    ):
        new_instance = cls(
            keep_sorted=keep_sorted,
            is_checked=is_checked,
            gene_id_to_gene_index={},
            transcript_ids_to_gene_ids_index={}
        )
        for feature in feature_iterator:
            new_instance = new_instance._add(feature)
        return new_instance

    def to_feature_iterator(self) -> Iterator[Feature]:
        for gene in self._gene_id_to_gene_index.values():
            yield gene
            for transcript in gene.transcript_values:
                yield transcript
                yield from transcript.exons

    @classmethod
    def from_feature_iterator(
            cls,
            feature_iterator: Iterable[Feature],
            keep_sorted: bool = False,
            is_checked: bool = False,
    ):
        feature_list = list(feature_iterator)
        exons = list(
            Exon(data=feature, is_checked=is_checked, shortcut=False) for feature in
            filter(lambda feature: feature.parsed_feature == FeatureType.EXON, feature_list)
        )
        transcripts = list(
            Transcript(
                data=feature,
                exons=[],
                is_checked=is_checked,
                is_inferred=False,
                keep_sorted=keep_sorted,
                shortcut=False
            ) for feature in
            filter(lambda feature: feature.parsed_feature == FeatureType.TRANSCRIPT, feature_list)
        )
        genes = list(
            Gene(
                data=feature,
                transcripts=[],
                transcript_ids=[],
                is_checked=is_checked,
                is_inferred=False,
                keep_sorted=keep_sorted,
                shortcut=False
            ) for feature in
            filter(lambda feature: feature.parsed_feature == FeatureType.GENE, feature_list)
        )
        transcript_ids = set(transcript.transcript_id for transcript in transcripts)
        for exon in tqdm(exons, desc="Scanning for missing transcript definitions..."):
            if exon.transcript_id not in transcript_ids:
                _lh.warning("Transcript %s inferred from exon!", exon.transcript_id)
                new_transcript = Transcript(
                    data=exon.get_data(),
                    exons=[],
                    is_checked=is_checked,
                    is_inferred=True,
                    keep_sorted=keep_sorted,
                    shortcut=False
                )
                transcripts.append(new_transcript)
                transcript_ids.add(new_transcript.transcript_id)
        transcript_id_to_transcript_index: Dict[str, Transcript] = {
            transcript.transcript_id: transcript
            for transcript in transcripts
        }
        for exon in tqdm(exons, desc="Adding exons to transcript..."):
            transcript_id_to_transcript_index[exon.transcript_id] = \
                transcript_id_to_transcript_index[exon.transcript_id].add_exon(exon)

        transcripts = transcript_id_to_transcript_index.values()
        gene_ids = set(gene.gene_id for gene in genes)
        for transcript in tqdm(transcripts, desc="Scanning for missing gene definitions..."):
            if transcript.gene_id not in gene_ids:
                _lh.warning("Gene %s inferred from transcript %s!", transcript.gene_id, transcript.transcript_id)
                new_gene = Gene(
                    data=transcript.get_data(),
                    transcripts=[],
                    transcript_ids=[],
                    is_checked=is_checked,
                    is_inferred=True,
                    keep_sorted=keep_sorted,
                    shortcut=False
                )
                genes.append(new_gene)
                gene_ids.add(new_gene.gene_id)
        gene_id_to_gene_index: Dict[str, Gene] = {
            gene.gene_id: gene
            for gene in genes
        }
        for transcript in tqdm(transcripts, desc="Adding transcripts to gene..."):
            gene_id_to_gene_index[transcript.gene_id] = \
                gene_id_to_gene_index[transcript.gene_id].add_transcript(transcript)
        return cls(
            gene_id_to_gene_index=gene_id_to_gene_index,
            transcript_ids_to_gene_ids_index={
                transcript.transcript_id: transcript.gene_id
                for transcript in transcripts
            },
            is_checked=is_checked,
            keep_sorted=keep_sorted
        )

    @classmethod
    def from_gtf_file(
            cls,
            gtf_file_path: str,
            keep_sorted: bool = False,
            is_checked: bool=False
    ):
        gtf_index_file_path = f"{gtf_file_path}.{GVPKL_VERSION}.gvpkl.xz"
        if should_regenerate(gtf_file_path, gtf_index_file_path):
            new_instance = cls.from_feature_iterator(
                GtfIterator(gtf_file_path),
                keep_sorted=keep_sorted,
                is_checked=is_checked
            )
            pickle_helper.dump((GVPKL_VERSION, new_instance), gtf_index_file_path)
            return new_instance
        else:
            (gvpkl_version, new_instance) = pickle_helper.load(gtf_index_file_path)
            if gvpkl_version != GVPKL_VERSION:
                new_instance = cls.from_feature_iterator(
                    GtfIterator(gtf_file_path),
                    keep_sorted=keep_sorted,
                    is_checked=is_checked
                )
                pickle_helper.dump((GVPKL_VERSION, new_instance), gtf_index_file_path)
            return new_instance


