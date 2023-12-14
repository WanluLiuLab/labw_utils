from __future__ import annotations

import gc
import os
import re
import tempfile

import msgpack
import pandas as pd

from blast_utils.blast6_io import parse_blast6, read_blast6
from blast_utils.blast_db import BlastSearchConfig, BlastConfig, AbstractBlastDB
from labw_utils.bioutils.datastructure.fai_view import create_fai_from_fasta
from labw_utils.bioutils.parser.fasta import FastaIterator
from labw_utils.commonutils.lwio.safe_io import get_reader, get_appender, get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import easyexec
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
from labw_utils.typing_importer import Dict, Mapping, Optional, Literal

_lh = get_logger(__name__)


class SingleBlastDB(AbstractBlastDB):
    _path: str
    _dbid: int
    _bc: BlastConfig
    _title_dict: Optional[Dict[str, str]]
    _txid_dict: Optional[Dict[str, int]]

    @property
    def bc(self) -> BlastConfig:
        return self._bc

    @property
    def dbid(self) -> int:
        return self._dbid

    @property
    def title_dict(self) -> Mapping[str, str]:
        if self._title_dict is None:
            with get_reader(os.path.join(self._path, "title.msgpack"), is_binary=True) as r:
                self._title_dict = msgpack.load(r)
        return self._title_dict

    @property
    def txid_dict(self) -> Mapping[str, int]:
        if self._txid_dict is None:
            with get_reader(os.path.join(self._path, "txid.msgpack"), is_binary=True) as r:
                self._txid_dict = msgpack.load(r)
        return self._txid_dict

    def __init__(
        self,
        tax_split_path: str,
        dbid: int,
        db_type: Literal["prot", "nucl"],
        bc: BlastConfig,
    ):
        self._dbid = dbid
        self._path = os.path.join(tax_split_path, str(self._dbid))
        self._db_type = db_type
        os.makedirs(self._path, exist_ok=True)
        self._bc = bc
        self._title_dict = None
        self._txid_dict = None

    def add(self, accession: str, txid: int, title: str, seq: str):
        safe_accession = re.sub(r"\W", "_", accession)
        safe_title = re.sub(r"\W", "_", title)
        with get_appender(os.path.join(self._path, "seq.fa")) as faa:
            faa.write(f">{safe_accession}\n{seq}\n")
        with get_appender(os.path.join(self._path, "title.tsv")) as tta:
            tta.write(f"{safe_accession}\t{safe_title}\n")
        with get_appender(os.path.join(self._path, "txid.tsv")) as tta:
            tta.write(f"{safe_accession}\t{txid}\n")

    def finalize(self):
        def _finalize_title_msgpack():
            self._title_dict = {}
            with get_reader(os.path.join(self._path, "title.tsv"), is_binary=False) as r:
                for l in r:
                    l = l.strip()
                    k, v = l.split("\t")
                    self._title_dict[k] = v
            with get_writer(os.path.join(self._path, "title.msgpack"), is_binary=True) as w:
                msgpack.dump(self._title_dict, w)
            rm_rf(os.path.join(self._path, "title.tsv"))

        def _finalize_txid_msgpack():
            self._txid_dict = {}
            with get_reader(os.path.join(self._path, "txid.tsv"), is_binary=False) as r:
                for l in r:
                    l = l.strip()
                    k, v = l.split("\t")
                    self._txid_dict[k] = int(v)
            with get_writer(os.path.join(self._path, "txid.msgpack"), is_binary=True) as w:
                msgpack.dump(self._txid_dict, w)
            rm_rf(os.path.join(self._path, "txid.tsv"))

        def _finalize_fai():
            if self._bc.samtools_path is not None:
                easyexec(
                    [
                        self._bc.samtools_path,
                        "faidx",
                        os.path.join(self._path, "seq.fa"),
                    ],
                    log_path=None,
                )
            else:
                _lh.warning("samtools not installed; Using slow implementation of FASTA indexer.")
                create_fai_from_fasta(
                    fasta_filename=os.path.join(self._path, "seq.fa"),
                    fai_filename=os.path.join(self._path, "seq.fa.fai"),
                )

        def _finalize_mmseqs():
            if self._bc.mmseqs_path is not None:
                easyexec(
                    [
                        self._bc.mmseqs_path,
                        "createdb",
                        os.path.join(self._path, "seq.fa"),
                        os.path.join(self._path, "mmseqs_db"),
                    ],
                    log_path=os.path.join(self._path, "mmseqs.log"),
                )
            else:
                _lh.warning("mmseqs not installed; Creation of mmseqs database skipped.")

        def _finalize_blast():
            blast_txid_map_path = os.path.join(self._path, "blast_txid_map.tsv")
            with get_writer(blast_txid_map_path, is_binary=False) as writer:
                for accession, txid in self.txid_dict.items():
                    writer.write(f"{accession}\t{txid}\n")

            if self._bc.makeblastdb_path is not None:
                easyexec(
                    [
                        self._bc.makeblastdb_path,
                        "-dbtype",
                        self._db_type,
                        "-blastdb_version",
                        "5",
                        "-in",
                        os.path.join(self._path, "seq.fa"),
                        "-out",
                        os.path.join(self._path, ""),
                        "-taxid_map",
                        blast_txid_map_path,
                        "-parse_seqids",
                    ],
                    log_path=os.path.join(self._path, "blast_db.log"),
                )
            else:
                _lh.warning("makeblastdb not installed; Creation of blast database skipped.")

        def _finalize_diamond():
            diamond_txid_map_path = os.path.join(self._path, "diamond_taxon.tsv")
            with get_writer(diamond_txid_map_path, is_binary=False) as writer:
                writer.write(f"accession.version\ttaxid\n")
                for accession, txid in self.txid_dict.items():
                    writer.write(f"{accession}\t{txid}\n")
            if self._bc.diamond_path is not None:
                cmd = [
                    self._bc.diamond_path,
                    "makedb",
                    "--threads",
                    "1",
                    "--in",
                    os.path.join(self._path, "seq.fa"),
                    "--db",
                    os.path.join(self._path, "diamond_db"),
                    "--ignore-warnings",  # Force build on truncated databases
                ]
                if self._bc.ncbi_txdb_names_path is not None and self._bc.ncbi_txdb_nodes_path is not None:
                    cmd.extend(
                        [
                            "--taxonmap",
                            diamond_txid_map_path,
                            "--taxonnodes",
                            self._bc.ncbi_txdb_nodes_path,
                            "--taxonnames",
                            self._bc.ncbi_txdb_names_path,
                        ]
                    )
                easyexec(
                    cmd,
                    log_path=os.path.join(self._path, "diamond_db.log"),
                )
            else:
                _lh.warning("diamond not installed; Creation of diamond database skipped.")

        if not os.path.exists(os.path.join(self._path, "title.msgpack")):
            _finalize_title_msgpack()
        if not os.path.exists(os.path.join(self._path, "txid.msgpack")):
            _finalize_txid_msgpack()
        if not os.path.exists(os.path.join(self._path, "seq.fa.fai")):
            _finalize_fai()
        if not os.path.exists(os.path.join(self._path, "mmseqs_db")):
            _finalize_mmseqs()
        if self._db_type == "prot" and not os.path.exists(os.path.join(self._path, "diamond_db.dmnd")):
            _finalize_diamond()
        if not os.path.exists(os.path.join(self._path, "blast_db.pjs")):
            _finalize_blast()
        self.gc()

    def search_exec(
        self,
        bsc: BlastSearchConfig,
        query_path: str,
        dest_path: str,
        tmpdir: Optional[str] = None,
    ):
        _ = tmpdir
        del tmpdir
        log_path = os.path.join(dest_path, f"{bsc.search_executable}_search.log")

        def _get_seqdict(
            blast6: pd.DataFrame,
        ) -> Mapping[str, str]:
            if len(blast6) == 0:
                return {}
            blast6_bed: pd.DataFrame = blast6.loc[:, ["sseqid", "sstart", "send"]]

            blast6_bed["sstart"] = blast6_bed["sstart"] - 1
            seq_dict = {}

            with tempfile.TemporaryDirectory() as tmpdirname:
                tmpbed = os.path.join(tmpdirname, "tmp.bed")
                tmpfa = os.path.join(tmpdirname, "tmp.fa")
                blast6_bed.to_csv(tmpbed, sep="\t", index=False, header=False)
                easyexec(
                    [
                        self._bc.bedtools_path,
                        "getfasta",
                        "-fi",
                        os.path.join(self._path, "seq.fa"),
                        "-bed",
                        tmpbed,
                        "-fo",
                        tmpfa,
                    ],
                    None,
                )
                with FastaIterator(tmpfa, show_tqdm=False) as fai:
                    for record in fai:
                        seq_dict[record.seq_id] = record.sequence
            return seq_dict

        def _search_exec_blast():
            if bsc.search_executable in {
                "blastn",
                "tblastn",
                "tblastx",
            }:
                if self._db_type != "nt":
                    raise RuntimeError(
                        f"BLAST program {bsc.search_executable} is designed for nucleotide, "
                        "but this BLAST database is built for protein"
                    )
            else:
                if self._db_type != "prot":
                    raise RuntimeError(
                        f"BLAST program {bsc.search_executable} is designed for protein, "
                        "but this BLAST database is built for nucleotide"
                    )
            easyexec(
                [
                    self._bc.which(bsc.search_executable),
                    "-query",
                    query_path,
                    "-db",
                    os.path.join(self._path, "blast_db"),
                    "-out",
                    blast_result_file,
                    "-outfmt",
                    str(7),
                    "-num_threads",
                    str(bsc.n_jobs),
                    "-evalue",
                    str(bsc.filter_e_value),
                    *bsc.cmdargs,
                ],
                log_path=log_path,
            )

        def _search_exec_mmseqs():
            with tempfile.TemporaryDirectory() as tmpdir:
                easyexec(
                    [
                        self._bc.mmseqs_path,
                        "easy-search",
                        query_path,
                        os.path.join(self._path, "mmseqs_db"),
                        blast_result_file,
                        tmpdir,
                        "--format-mode",
                        str(0),
                        "--threads",
                        str(bsc.n_jobs),
                        *bsc.cmdargs,
                    ],
                    log_path=log_path,
                )

        def _search_exec_diamond():
            easyexec(
                [
                    *bsc.search_executable.split("_"),
                    "--query",
                    query_path,
                    "--db",
                    os.path.join(self._path, "diamond_db"),
                    "--out",
                    blast_result_file,
                    "--outfmt",
                    str(6),
                    "--evalue",
                    str(bsc.filter_e_value),
                    "--threads",
                    str(bsc.n_jobs),
                    *bsc.cmdargs,
                ],
                log_path=log_path,
            )

        blast_result_file = os.path.join(dest_path, "search_result.blast6.tsv")
        if bsc.search_executable.startswith("diamond"):
            self._bc.require(["diamond"])
        else:
            self._bc.require([bsc.search_executable])
        if bsc.search_executable in {
            "blastp",
            "blastx",
            "blastn",
            "psiblast",
            "tblastn",
            "tblastx",
        }:
            _search_exec_blast()
        elif bsc.search_executable == "mmseqs":
            _search_exec_mmseqs()
        elif bsc.search_executable in {"diamond_blastp", "diamond_blastx"}:
            if self._db_type != "prot":
                raise RuntimeError(
                    f"BLAST program {bsc.search_executable} is designed for protein, "
                    "but this BLAST database is built for nucleotide"
                )
            _search_exec_diamond()
        else:
            raise ValueError(f"Unknown program {bsc.search_executable}")
        result_blast6 = read_blast6(blast_result_file)
        parse_blast6(
            src_blast6_df=result_blast6,
            title_dict=self.title_dict,
            txid_dict=self.txid_dict,
            seq_dict=_get_seqdict(result_blast6),
            bsc=bsc,
        ).to_csv(
            os.path.join(dest_path, "search_result.parsed.tsv"),
            sep="\t",
            index=False,
        )
        self.gc()

    def gc(self):
        self._title_dict = None
        self._txid_dict = None
        gc.collect()
