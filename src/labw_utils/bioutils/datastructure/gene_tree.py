"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

import itertools
from abc import abstractmethod, ABC
from collections import defaultdict

from labw_utils.bioutils.datastructure.feature_view import FeatureView
from labw_utils.bioutils.datastructure.gv import GVPError, CanCheckInterface
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.datastructure.gv.gene import Gene, DumbGene
from labw_utils.bioutils.datastructure.gv.gene_container_interface import (
    GeneContainerInterface,
)
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container_interface import (
    TranscriptContainerInterface,
    DuplicatedTranscriptIDError,
)
from labw_utils.bioutils.record.feature import Feature, FeatureInterface, FeatureType
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.file_system import should_regenerate
from labw_utils.commonutils.stdlib_helper import pickle_helper
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import (
    Iterable,
    Dict,
    Iterator,
    Sequence,
    Optional,
    Mapping,
    Literal,
    Sized,
)

from labw_utils.typing_importer import SequenceProxy


_lh = get_logger(__name__)

GVPKL_VERSION = "1.2.1"
"""
Current version of GVPKL standard.

Changes:

- 1.1: Added gc() and multiple transcribe_* alternative to Transcript.
- 1.2: Added BiologicalInterval class.
- 1.2.1: Bugs in collapsing fixed

.. versionadded:: 1.0.3
"""


class DuplicatedGeneIDError(GVPError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    def __init__(self, gene_id: str):
        super().__init__(f"Gene ID {gene_id} duplicated")


class GeneTreeInterface(GeneContainerInterface, TranscriptContainerInterface, CanCheckInterface, Sized):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    @abstractmethod
    def add_exon(self, exon: Exon) -> GeneTreeInterface:
        raise NotImplementedError

    @abstractmethod
    def del_exon(self, transcript_id: str, exon_index: int) -> GeneTreeInterface:
        raise NotImplementedError

    @abstractmethod
    def to_feature_iterator(self) -> Iterator[FeatureInterface]:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_feature_iterator(
        cls,
        feature_iterator: Iterable[Feature],
        shortcut: bool = False,
        is_checked: bool = False,
        show_tqdm: bool = True,
        gene_implementation: Literal[Gene, DumbGene] = Gene,
    ) -> GeneTreeInterface:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_gtf_file(cls, gtf_file_path: str, show_tqdm: bool = True, is_checked: bool = False) -> GeneTreeInterface:
        raise NotImplementedError

    @classmethod
    def from_gvpkl(cls, gtf_index_file_path: str) -> GeneTreeInterface:
        raise NotImplementedError

    def to_gvpkl(self, gtf_index_file_path: str):
        raise NotImplementedError

    @abstractmethod
    def gc(self):
        raise NotImplementedError


class BaseGeneTree(GeneTreeInterface, ABC):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    @classmethod
    def from_gvpkl(cls, gtf_index_file_path: str):
        (gvpkl_version, new_instance) = pickle_helper.load(gtf_index_file_path)
        if gvpkl_version == GVPKL_VERSION:
            return new_instance
        else:
            raise TypeError(f"Version {gvpkl_version} (file) and {GVPKL_VERSION} (library) mismatch")

    def to_gvpkl(self, gtf_index_file_path: str):
        pickle_helper.dump((GVPKL_VERSION, self), gtf_index_file_path)

    @classmethod
    def from_gtf_file(
        cls,
        gtf_file_path: str,
        is_checked: bool = False,
        show_tqdm: bool = True,
        gene_implementation: Literal[Gene, DumbGene] = Gene,
    ) -> GeneTreeInterface:
        gtf_index_file_path = f"{gtf_file_path}.{GVPKL_VERSION}.gvpkl.xz"
        if not should_regenerate(gtf_file_path, gtf_index_file_path):
            try:
                return cls.from_gvpkl(gtf_index_file_path)
            except (TypeError, EOFError):
                pass
        new_instance = cls.from_feature_iterator(
            FeatureView.from_gtf(gtf_file_path, show_tqdm=show_tqdm),
            is_checked=is_checked,
            gene_implementation=gene_implementation,
            show_tqdm=show_tqdm,
        )
        new_instance.to_gvpkl(gtf_index_file_path)
        return new_instance

    def gc(self):
        for gene in self.gene_values:
            gene.gc()


class GeneTree(BaseGeneTree):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    __slots__ = (
        "_gene_id_to_gene_index",
        "_transcript_ids_to_gene_ids_index",
        "_transcripts",
    )
    _gene_id_to_gene_index: Dict[str, Gene]
    _transcript_ids_to_gene_ids_index: Dict[str, str]
    """transcript_id -> gene_id"""

    def __len__(self) -> int:
        i = 0
        for gene in self._gene_id_to_gene_index.values():
            i += 1
            for transcript in gene.transcript_values:
                i += 1
                i += transcript.number_of_exons
        return i

    _transcripts: Optional[Sequence[Transcript]]

    @property
    def number_of_genes(self) -> int:
        return len(self._gene_id_to_gene_index)

    @property
    def gene_values(self) -> Sequence[Gene]:
        return SequenceProxy(self._gene_id_to_gene_index.values())

    @property
    def gene_ids(self) -> Sequence[str]:
        return SequenceProxy(self._gene_id_to_gene_index.keys())

    @property
    def number_of_transcripts(self) -> int:
        return len(self.transcript_values)

    @property
    def transcript_values(self) -> Sequence[Transcript]:
        if self._transcripts is None:
            self._transcripts = tuple(self.get_transcript(transcript_id) for transcript_id in self.transcript_ids)
        return SequenceProxy(self._transcripts)

    @property
    def transcript_ids(self) -> Sequence[str]:
        return SequenceProxy(self._transcript_ids_to_gene_ids_index.keys())

    def __init__(
        self,
        *,
        is_checked: bool,
        shortcut: bool,
        gene_id_to_gene_index: Mapping[str, Gene],
        transcript_ids_to_gene_ids_index: Mapping[str, str],
    ):
        self._is_checked = is_checked
        if not shortcut:
            self._gene_id_to_gene_index = dict(gene_id_to_gene_index)
            self._transcript_ids_to_gene_ids_index = dict(transcript_ids_to_gene_ids_index)
        else:
            self._gene_id_to_gene_index = gene_id_to_gene_index  # type: ignore
            self._transcript_ids_to_gene_ids_index = transcript_ids_to_gene_ids_index  # type: ignore
        self._transcripts = None

    def get_gene(self, gene_id: str) -> Sequence[Gene]:
        try:
            return SequenceProxy([self._gene_id_to_gene_index[gene_id]])
        except KeyError:
            return SequenceProxy.empty()

    def add_gene(self, gene: Gene) -> GeneTree:
        if not self._is_checked and gene.gene_id in self._gene_id_to_gene_index:
            raise DuplicatedGeneIDError(gene.gene_id)
        new_gene_id_to_gene_index = dict(self._gene_id_to_gene_index)
        new_transcript_ids_to_gene_ids_index = dict(self._transcript_ids_to_gene_ids_index)
        new_gene_id_to_gene_index[gene.gene_id] = gene
        for transcript_id in gene.transcript_ids:
            new_transcript_ids_to_gene_ids_index[transcript_id] = gene.gene_id
        return GeneTree(
            is_checked=self._is_checked,
            shortcut=True,
            gene_id_to_gene_index=new_gene_id_to_gene_index,
            transcript_ids_to_gene_ids_index=new_transcript_ids_to_gene_ids_index,
        )

    def del_gene(self, gene_id: str) -> GeneTree:
        new_gene_id_to_gene_index = dict(self._gene_id_to_gene_index)
        _ = new_gene_id_to_gene_index.pop(gene_id)
        new_transcript_ids_to_gene_ids_index = {
            k: v for k, v in self._transcript_ids_to_gene_ids_index.items() if v != gene_id
        }
        return GeneTree(
            is_checked=self._is_checked,
            shortcut=True,
            gene_id_to_gene_index=new_gene_id_to_gene_index,
            transcript_ids_to_gene_ids_index=new_transcript_ids_to_gene_ids_index,
        )

    def replace_gene(self, new_gene: Gene) -> GeneTree:
        return self.del_gene(new_gene.gene_id).add_gene(new_gene)

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self.get_gene(self._transcript_ids_to_gene_ids_index[transcript_id])[0].get_transcript(transcript_id)

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
                    shortcut=False,
                    transcripts=[],
                    transcript_ids=[],
                )
            ).add_transcript(transcript)

    def del_transcript(self, transcript_id: str) -> GeneTree:
        gene_to_be_modified = self._gene_id_to_gene_index[self._transcript_ids_to_gene_ids_index[transcript_id]]
        gene_to_be_modified = gene_to_be_modified.del_transcript(transcript_id)
        self._gene_id_to_gene_index[gene_to_be_modified.gene_id] = gene_to_be_modified
        _ = self._transcript_ids_to_gene_ids_index.pop(transcript_id)
        return self

    def replace_transcript(self, new_transcript: Transcript) -> GeneTree:
        return self.del_transcript(new_transcript.transcript_id).add_transcript(new_transcript)

    def add_exon(self, exon: Exon) -> GeneTree:
        if exon.transcript_id in self._transcript_ids_to_gene_ids_index:
            return self.replace_transcript(self.get_transcript(exon.transcript_id).add_exon(exon))
        else:
            return self.add_transcript(
                Transcript(
                    data=exon.get_data(),
                    exons=[],
                    is_inferred=True,
                    is_checked=self._is_checked,
                    shortcut=False,
                )
            ).add_exon(exon)

    def del_exon(self, transcript_id: str, exon_index: int) -> GeneTree:
        return self.replace_transcript(self.get_transcript(transcript_id).del_exon(exon_index))

    def to_feature_iterator(self) -> Iterator[FeatureInterface]:
        for gene in self._gene_id_to_gene_index.values():
            if not isinstance(gene, DumbGene):
                yield gene
            for transcript in gene.transcript_values:
                yield transcript
                yield from transcript.exons

    @classmethod
    def from_feature_iterator(
        cls,
        feature_iterator: Iterable[FeatureInterface],
        shortcut: bool = False,
        is_checked: bool = False,
        show_tqdm: bool = True,
        gene_implementation: Literal[Gene, DumbGene] = Gene,
    ):
        if not shortcut:
            feature_list = list(feature_iterator)
        else:
            feature_list = feature_iterator
        _lh.info("Filtering for gene, transcript and exon definition...")
        initially_added_exons = list(
            Exon(data=feature, is_checked=is_checked, shortcut=False)
            for feature in filter(lambda feature: feature.parsed_feature == FeatureType.EXON, feature_list)
        )
        initially_added_transcripts = list(
            Transcript(
                data=feature,
                exons=[],
                is_checked=is_checked,
                is_inferred=False,
                shortcut=False,
            )
            for feature in filter(
                lambda feature: feature.parsed_feature == FeatureType.TRANSCRIPT,
                feature_list,
            )
        )
        initially_added_genes = list(
            gene_implementation(
                data=feature,
                transcripts=[],
                transcript_ids=[],
                is_checked=is_checked,
                is_inferred=False,
                shortcut=False,
            )
            for feature in filter(lambda feature: feature.parsed_feature == FeatureType.GENE, feature_list)
        )

        # Normalize transcript IDs

        transcript_ids = set()
        it = initially_added_transcripts
        if show_tqdm:
            it = tqdm(it, desc="Scanning for duplicated transcript definitions...")
        for transcript in it:
            if transcript.transcript_id in transcript_ids:
                raise DuplicatedTranscriptIDError(transcript.transcript_id)
            transcript_ids.add(transcript.transcript_id)
        it = initially_added_exons
        if show_tqdm:
            it = tqdm(it, desc="Scanning for missing transcript definitions...")
        for exon in it:
            if exon.transcript_id not in transcript_ids:
                _lh.warning("Transcript %s inferred from exon!", exon.transcript_id)
                new_transcript = Transcript(
                    data=exon.get_data(),
                    exons=[],
                    is_checked=is_checked,
                    is_inferred=True,
                    shortcut=False,
                )
                initially_added_transcripts.append(new_transcript)
                transcript_ids.add(new_transcript.transcript_id)

        transcript_id_to_transcript_index: Dict[str, Transcript] = {
            transcript.transcript_id: transcript for transcript in initially_added_transcripts
        }
        if show_tqdm:
            it = tqdm(it, desc="Adding exons to transcript...")
        for exon in it:
            transcript_id_to_transcript_index[exon.transcript_id] = transcript_id_to_transcript_index[
                exon.transcript_id
            ].add_exon(exon)

        finalized_transcripts = transcript_id_to_transcript_index.values()

        # Normalize Gene IDs

        gene_ids = set()
        it = initially_added_genes
        if show_tqdm:
            it = tqdm(it, desc="Scanning for duplicated gene definitions...")
        for gene in it:
            if gene.gene_id in gene_ids:
                raise DuplicatedGeneIDError(gene.gene_id)
            gene_ids.add(gene.gene_id)
        it = finalized_transcripts
        if show_tqdm:
            it = tqdm(it, desc="Scanning for missing gene definitions...")
        for transcript in it:
            if transcript.gene_id not in gene_ids:
                # _lh.warning("Gene %s inferred from transcript %s!", transcript.gene_id, transcript.transcript_id)
                new_gene = gene_implementation(
                    data=transcript.get_data(),
                    transcripts=[],
                    transcript_ids=[],
                    is_checked=is_checked,
                    is_inferred=True,
                    shortcut=False,
                )
                initially_added_genes.append(new_gene)
                gene_ids.add(new_gene.gene_id)
        gene_id_to_gene_index: Dict[str, Gene] = {gene.gene_id: gene for gene in initially_added_genes}
        it = finalized_transcripts
        if show_tqdm:
            it = tqdm(it, desc="Scanning for missing gene definitions...")
        for transcript in it:
            gene_id_to_gene_index[transcript.gene_id] = gene_id_to_gene_index[transcript.gene_id].add_transcript(
                transcript
            )

        return cls(
            gene_id_to_gene_index=gene_id_to_gene_index,
            transcript_ids_to_gene_ids_index={
                transcript.transcript_id: transcript.gene_id for transcript in finalized_transcripts
            },
            shortcut=True,
            is_checked=is_checked,
        )


class DiploidGeneTree(BaseGeneTree):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    def __len__(self) -> int:
        return sum(len(gt) for gt in self._chr_gt_idx.values())

    _chr_gt_idx: Dict[str, GeneTreeInterface]

    def __init__(
        self,
        *,
        is_checked: bool,
        shortcut: bool,
        chr_gt_idx: Mapping[str, GeneTreeInterface],
    ):
        self._is_checked = is_checked
        if not shortcut:
            self._chr_gt_idx = dict(chr_gt_idx)
        else:
            self._chr_gt_idx = chr_gt_idx  # type: ignore

    def _get_or_create_gt(self, feature: BaseFeatureProxy):
        if feature.seqname not in self._chr_gt_idx:
            self._chr_gt_idx[feature.seqname] = GeneTree(
                gene_id_to_gene_index={},
                transcript_ids_to_gene_ids_index={},
                shortcut=True,
                is_checked=self.is_checked,
            )
        return self._chr_gt_idx[feature.seqname]

    def add_exon(self, exon: Exon) -> GeneTreeInterface:
        this_gt = self._get_or_create_gt(exon)
        this_gt = this_gt.add_exon(exon)
        self._chr_gt_idx[exon.seqname] = this_gt
        return self

    def del_exon(self, transcript_id: str, exon_index: int) -> GeneTreeInterface:
        raise NotImplementedError

    def to_feature_iterator(self) -> Iterator[FeatureInterface]:
        for gt in self._chr_gt_idx.values():
            yield from gt.to_feature_iterator()

    @classmethod
    def from_feature_iterator(
        cls,
        feature_iterator: Iterable[FeatureInterface],
        shortcut: bool = False,
        is_checked: bool = False,
        show_tqdm: bool = True,
        gene_implementation: Literal[Gene, DumbGene] = Gene,
    ):
        _ = shortcut
        del shortcut
        feature_chr_dict = defaultdict(lambda: [])
        for feature in feature_iterator:
            feature_chr_dict[feature.seqname].append(feature)
        chr_gt_idx = {}
        it = feature_chr_dict.items()
        if show_tqdm:
            it = tqdm(it, desc="Processing contigs...")
        for k, v in it:
            chr_gt_idx[k] = GeneTree.from_feature_iterator(
                v,
                shortcut=True,
                is_checked=is_checked,
                gene_implementation=gene_implementation,
                show_tqdm=False,
            )
        return cls(chr_gt_idx=chr_gt_idx, is_checked=is_checked, shortcut=True)

    @property
    def number_of_genes(self) -> int:
        return sum(gt.number_of_genes for gt in self._chr_gt_idx.values())

    @property
    def gene_values(self) -> Sequence[Gene]:
        return SequenceProxy(itertools.chain(itertools.chain(*(gt.gene_values for gt in self._chr_gt_idx.values()))))

    @property
    def gene_ids(self) -> Sequence[str]:
        return SequenceProxy(list(set(itertools.chain(*(gt.gene_ids for gt in self._chr_gt_idx.values())))))

    def get_gene(self, gene_id: str) -> Sequence[Gene]:
        return SequenceProxy(list(itertools.chain(*(gt.get_gene(gene_id) for gt in self._chr_gt_idx.values()))))

    def add_gene(self, gene: Gene) -> GeneContainerInterface:
        this_gt = self._get_or_create_gt(gene)
        this_gt = this_gt.add_gene(gene)
        self._chr_gt_idx[gene.seqname] = this_gt
        return self

    def del_gene(self, gene_id: str) -> GeneContainerInterface:
        raise NotImplementedError

    def replace_gene(self, new_gene: Gene) -> GeneContainerInterface:
        raise NotImplementedError

    @property
    def number_of_transcripts(self) -> int:
        return sum(gt.number_of_transcripts for gt in self._chr_gt_idx.values())

    @property
    def transcript_values(self) -> Sequence[Transcript]:
        return SequenceProxy(list(itertools.chain(*(gt.transcript_values for gt in self._chr_gt_idx.values()))))

    @property
    def transcript_ids(self) -> Sequence[str]:
        return SequenceProxy(list(itertools.chain(*(gt.transcript_ids for gt in self._chr_gt_idx.values()))))

    def get_transcript(self, transcript_id: str) -> Transcript:
        raise NotImplementedError

    def add_transcript(self, transcript: Transcript) -> TranscriptContainerInterface:
        this_gt = self._get_or_create_gt(transcript)
        this_gt = this_gt.add_transcript(transcript)
        self._chr_gt_idx[transcript.seqname] = this_gt
        return self

    def del_transcript(self, transcript_id: str) -> TranscriptContainerInterface:
        raise NotImplementedError

    def replace_transcript(self, new_transcript: Transcript) -> TranscriptContainerInterface:
        raise NotImplementedError
