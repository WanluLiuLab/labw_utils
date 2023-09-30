from __future__ import annotations

import json
import math
import os.path
import re
import shutil
from typing import Optional, Sequence

import networkx as nx
import pandas as pd
from networkx.readwrite import json_graph

from blast_utils.blast_db import BlastSearchConfig, AbstractBlastDB
from blast_utils.blast_db._sbdb import SingleBlastDB
from blast_utils.ncbi_taxdb import NON_EXIST
from labw_utils.bioutils.parser.fasta import FastaIterator, FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.lwio import get_writer, get_reader
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
from labw_utils.typing_importer import Tuple, Set, Mapping, Any, List, Dict

_lh = get_logger(__name__)


def search_exec_rbh(
    dbdb: AbstractBlastDB,
    query_path: str,
    dest_path: str,
    bsc: BlastSearchConfig,
    tmpdir: Optional[str] = None,
):
    if not bsc.keep_longest_hit_only:
        raise ValueError("Please set `keep_longest_hit_only` to True.")
    _lh.info("Adding query sequence to current BLASTDB...")
    if tmpdir is None:
        tmpdir = os.path.join(dest_path, "tmp")
    rm_rf(tmpdir)  # Must be empty
    os.makedirs(tmpdir)

    # Create query BLAST DB
    query_blastdb = SingleBlastDB(
        os.path.join(tmpdir, "query_blastdb"),
        0,
        "prot",
        dbdb.bc,
    )
    with FastaIterator(query_path, show_tqdm=False) as fai:
        for record in fai:
            accession = re.sub(r"\W", "_", record.seq_id.split(" ")[0])
            # Query sequence should not have a title.
            query_blastdb.add(accession, NON_EXIST, "UNTITLED", record.sequence)
    query_blastdb.finalize()

    # Perform first-pass search
    first_pass_dest = os.path.join(tmpdir, "first_pass")
    first_pass_bsc = bsc.duplicate()
    dbdb.search_exec(
        bsc=first_pass_bsc,
        query_path=query_path,
        dest_path=first_pass_dest,
    )

    # Split results
    first_pass_results = pd.read_csv(
        os.path.join(first_pass_dest, "search_result.parsed.tsv"),
        sep="\t",
    )
    first_pass_query_accession: Set[Tuple[str, str]] = set(
        (t.qseqid, t.sseqid) for t in first_pass_results.itertuples(index=False)
    )
    first_pass_fasta_file_path = os.path.join(first_pass_dest, "match.fa")
    with FastaWriter(first_pass_fasta_file_path) as faw:
        for t in first_pass_results.itertuples(index=False):
            faw.write(
                # Encode query sequence IDs in case multiple query have same hit
                FastaRecord(seq_id="-".join((t.qseqid, t.sseqid)), sequence=t.seq)
            )
    _lh.info(
        "1st search finished with %d hits left",
        len(first_pass_query_accession),
    )

    # Perform second-pass search
    second_pass_dest = os.path.join(tmpdir, "second_pass")
    second_pass_bsc = bsc.duplicate()
    query_blastdb.search_exec(
        query_path=first_pass_fasta_file_path,
        bsc=second_pass_bsc,
        dest_path=second_pass_dest,
    )
    second_pass_results = pd.read_csv(
        os.path.join(second_pass_dest, "search_result.parsed.tsv"),
        sep="\t",
    )
    second_pass_query_accession = set(
        (t.sseqid, t.qseqid.split("-")[1]) for t in second_pass_results.itertuples(index=False)
    )
    final_query_accession = first_pass_query_accession.intersection(second_pass_query_accession)
    _lh.info(
        "2nd search finished with %d accessions left",
        len(final_query_accession),
    )
    final_query_indexes = [
        i
        for i, t in enumerate(first_pass_results.itertuples(index=False))
        if (t.qseqid, t.sseqid) in final_query_accession
    ]
    final_results: pd.DataFrame = first_pass_results.iloc[final_query_indexes, :]
    final_results.to_csv(
        os.path.join(dest_path, "search_result.parsed.tsv"),
        sep="\t",
        index=False,
    )
    dbdb.gc()


def search_exec_expansion(
    dbdb: AbstractBlastDB,
    query_path: str,
    dest_path: str,
    bsc: BlastSearchConfig,
    max_n_iters: Optional[int] = 10,
    max_n_seqs: Optional[int] = 50000,
    tmpdir: Optional[str] = None,
    resume: bool = False,
) -> Tuple[int, ExpansiveBlastResult]:
    if tmpdir is None:
        tmpdir = os.path.join(dest_path, "tmp")
    if not resume:
        rm_rf(tmpdir)  # Must be empty
    os.makedirs(tmpdir, exist_ok=True)
    if max_n_iters is None:
        max_n_iters = math.inf
    if max_n_seqs is None:
        max_n_seqs = math.inf
    n_iters = 1
    ebr = ExpansiveBlastResult.start_from_query_sequence(query_path)
    this_iter_ckp_path = None
    while n_iters < max_n_iters:
        this_iter_dst_path = os.path.join(tmpdir, f"pass_{n_iters}")
        this_iter_ckp_path = os.path.join(this_iter_dst_path, f"ckp.json")
        if resume and os.path.exists(this_iter_ckp_path):
            with get_reader(this_iter_ckp_path, is_binary=True) as r:
                ebr = ExpansiveBlastResult.from_dict(json.load(r))
                _lh.info(
                    "Skipped search iter %d",
                    n_iters,
                )
                if not ebr.seq_id_next_iter:
                    _lh.info("Converged")
                    break
                if len(ebr.ebrs) > max_n_seqs:
                    _lh.info("Number of sequence overflow")
                    break
                n_iters += 1
                continue

        _lh.info(
            "Start search iter %d with %d new sequences",
            n_iters,
            len(ebr.seq_id_next_iter),
        )
        this_iter_query_path = os.path.join(this_iter_dst_path, "match.fa")
        ebr.export_next_iter_to_fasta(this_iter_query_path)
        dbdb.search_exec(
            query_path=this_iter_query_path,
            dest_path=this_iter_dst_path,
            bsc=bsc.duplicate(),
            tmpdir=None,
        )
        this_iter_results = pd.read_csv(
            os.path.join(this_iter_dst_path, "search_result.parsed.tsv"),
            sep="\t",
        )
        is_converged = ebr.add_parsed_blast6(this_iter_results, added_in=n_iters)
        with get_writer(this_iter_ckp_path, is_binary=False) as w:
            json.dump(ebr.to_dict(), w, indent=4)
        if is_converged:
            _lh.info("Converged")
            break
        if len(ebr.ebrs) > max_n_seqs:
            _lh.info("Number of sequence overflow")
            break
        n_iters += 1
    if this_iter_ckp_path is not None:
        shutil.copy(this_iter_ckp_path, os.path.join(dest_path, "ckp.json"))
    return n_iters, ebr


class ExpansiveBlastRecord:
    ebr_id: int
    sequence: str
    added_in: int
    coordinates: Set[Tuple[str, int, int]]
    """[Accession, start, end]"""

    def __init__(
        self,
        ebr_id: int,
        sequence: str,
        added_in: int,
        coordinates: Optional[Set[Tuple[str, int, int]]] = None,
    ):
        if coordinates is None:
            coordinates = set()
        self.ebr_id = ebr_id
        self.sequence = sequence
        self.coordinates = coordinates
        self.added_in = added_in

    def to_dict(self):
        return {
            "ebr_id": self.ebr_id,
            "sequence": self.sequence,
            "coordinates": list(self.coordinates),
            "added_in": self.added_in,
        }

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(
            ebr_id=d["ebr_id"],
            sequence=d["sequence"],
            coordinates=set(map(tuple, d["coordinates"])),
            added_in=d["added_in"],
        )


class ExpansiveBlastResult:
    graph: nx.DiGraph
    ebrs: List[ExpansiveBlastRecord]
    seq_ebrid_map: Dict[str, int]
    blacklist_seqs: Set[str]
    accession_metadata: Dict[str, Tuple[str, int]]
    """Accession -> [Title, TXID]"""
    seq_id_next_iter: Set[int]

    def to_dict(self):
        return {
            "graph": json_graph.adjacency_data(self.graph),
            "ebrs": list(map(ExpansiveBlastRecord.to_dict, self.ebrs)),
            "seq_id_next_iter": list(self.seq_id_next_iter),
            "accession_metadata": self.accession_metadata,
            "seq_ebrid_map": self.seq_ebrid_map,
            "blacklist_seqs": list(self.blacklist_seqs),
        }

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(
            graph=json_graph.adjacency_graph(d["graph"], directed=True),
            ebrs=list(map(ExpansiveBlastRecord.from_dict, d["ebrs"])),
            seq_id_next_iter=d["seq_id_next_iter"],
            accession_metadata=d["accession_metadata"],
            seq_ebrid_map=d["seq_ebrid_map"],
            blacklist_seqs=d["blacklist_seqs"],
        )

    @classmethod
    def start_from_query_sequence(cls, src_fasta_file_path: str):
        seq_ebrid_map = {}
        graph = nx.DiGraph()
        ebrs: List[ExpansiveBlastRecord] = []
        accession_metadata: Dict[str, Tuple[str, int]] = {}
        seq_id_next_iter: Set[int] = set()
        with FastaIterator(src_fasta_file_path, show_tqdm=False) as fai:
            for record in fai:
                accession = re.sub(r"\W", "_", record.seq_id.split(" ")[0])
                seq = record.sequence
                accession_metadata[accession] = (accession, NON_EXIST)
                if seq not in seq_ebrid_map:
                    this_ebr_id = len(seq_ebrid_map)
                    seq_ebrid_map[seq] = this_ebr_id
                    ebrs.append(ExpansiveBlastRecord(this_ebr_id, seq, 0))
                    graph.add_node(this_ebr_id)
                    seq_id_next_iter.add(this_ebr_id)
                else:
                    this_ebr_id = seq_ebrid_map[seq]
                ebrs[this_ebr_id].coordinates.add((accession, 0, len(seq) - 1))

        return cls(
            graph=graph,
            ebrs=ebrs,
            accession_metadata=accession_metadata,
            seq_id_next_iter=seq_id_next_iter,
            seq_ebrid_map=seq_ebrid_map,
            blacklist_seqs=[],
        )

    def __init__(
        self,
        graph: nx.DiGraph,
        ebrs: Sequence[ExpansiveBlastRecord],
        seq_ebrid_map: Mapping[str, int],
        accession_metadata: Mapping[str, Tuple[str, int]],
        seq_id_next_iter: Set[int],
        blacklist_seqs: Sequence[str],
    ):
        self.graph = graph
        self.ebrs = list(ebrs)
        self.seq_ebrid_map = dict(seq_ebrid_map)
        self.accession_metadata = dict(accession_metadata)
        self.seq_id_next_iter = set(seq_id_next_iter)
        self.blacklist_seqs = set(blacklist_seqs)

    def export_next_iter_to_fasta(self, dst_fasta_file_path: str):
        with FastaWriter(dst_fasta_file_path) as faw:
            for i in self.seq_id_next_iter:
                faw.write(FastaRecord(seq_id=str(i), sequence=self.ebrs[i].sequence))

    def add_parsed_blast6(self, parsed_blast6: pd.DataFrame, added_in: int) -> bool:
        """
        :param parsed_blast6: Dataframe of parsed BLAST6
        :return: Whether the search had converged.
        """
        self.seq_id_next_iter = set()
        for tup in parsed_blast6.itertuples(index=False):
            dst_acc, dst_seq = tup.sseqid, tup.seq
            if dst_seq in self.blacklist_seqs:
                continue
            if dst_acc not in self.accession_metadata:
                self.accession_metadata[dst_acc] = (tup.title, tup.txid)
            if dst_seq not in self.seq_ebrid_map:
                dst_ebr_id = len(self.seq_ebrid_map)
                self.seq_ebrid_map[dst_seq] = dst_ebr_id
                self.ebrs.append(ExpansiveBlastRecord(dst_ebr_id, dst_seq, added_in=added_in))
                self.seq_id_next_iter.add(dst_ebr_id)
            else:
                dst_ebr_id = self.seq_ebrid_map[dst_seq]
            dst_start, dst_end = tup.sstart, tup.send
            self.ebrs[dst_ebr_id].coordinates.add((dst_acc, dst_start, dst_end))
            src_ebr_id = int(tup.qseqid)
            self.graph.add_edge(src_ebr_id, dst_ebr_id)

        return not self.seq_id_next_iter
