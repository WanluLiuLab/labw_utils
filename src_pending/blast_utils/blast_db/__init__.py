from __future__ import annotations

import abc
import multiprocessing
import shutil
from abc import abstractmethod

from blast_utils import AbstractConfig
from labw_utils.typing_importer import List, Literal, Optional

CHUNK_SIZE_DEFAULT = 1 << 14


class BlastSearchConfig(AbstractConfig):
    n_jobs: int
    search_executable: Literal[
        "blastp",
        "blastx",
        "blastn",
        "tblastn",
        "tblastx",
        "psiblast",
        "mmseqs",
        "diamond_blastp",
        "diamond_blastx",
    ]
    cmdargs: List[str]
    keep_longest_hit_only: bool
    filter_e_value: float
    remote_blast_db_name: str
    remote_blast_n_max_retries: int
    pident_cutoff: float

    def __init__(
        self,
        search_executable: Literal[
            "blastp",
            "blastx",
            "blastn",
            "tblastn",
            "tblastx",
            "psiblast",
            "mmseqs",
            "diamond_blastp",
            "diamond_blastx",
        ],
        n_jobs: int = multiprocessing.cpu_count(),
        cmdargs: Optional[List[str]] = None,
        keep_longest_hit_only: bool = False,
        filter_e_value: float = 1e-6,
        pident_cutoff: float = 80,
        remote_blast_db_name: str = "nr",
        remote_blast_n_max_retries: int = 10,
    ):
        self.n_jobs = n_jobs
        self.search_executable = search_executable
        self.cmdargs = cmdargs if cmdargs is not None else []
        self.keep_longest_hit_only = keep_longest_hit_only
        self.filter_e_value = filter_e_value
        self.remote_blast_db_name = remote_blast_db_name
        self.remote_blast_n_max_retries = remote_blast_n_max_retries
        self.pident_cutoff = pident_cutoff

    def duplicate(self) -> BlastSearchConfig:
        return self.__class__(
            search_executable=self.search_executable,
            n_jobs=self.n_jobs,
            cmdargs=self.cmdargs,
            keep_longest_hit_only=self.keep_longest_hit_only,
            remote_blast_db_name=self.remote_blast_db_name,
            remote_blast_n_max_retries=self.remote_blast_n_max_retries,
            pident_cutoff=self.pident_cutoff,
        )


class BlastConfig(AbstractConfig):
    mmseqs_path: Optional[str]
    diamond_path: Optional[str]
    blastn_path: Optional[str]
    blastp_path: Optional[str]
    blastx_path: Optional[str]
    tblastn_path: Optional[str]
    tblastx_path: Optional[str]
    psiblast_path: Optional[str]
    samtools_path: Optional[str]
    bedtools_path: Optional[str]
    makeblastdb_path: Optional[str]
    chunk_size: int
    ncbi_txdb_names_path: Optional[str]
    ncbi_txdb_nodes_path: Optional[str]

    def which(self, what: str) -> Optional[str]:
        return {
            "blastp": self.blastp_path,
            "blastx": self.blastx_path,
            "blastn": self.blastn_path,
            "tblastn": self.tblastn_path,
            "tblastx": self.tblastx_path,
            "mmseqs": self.mmseqs_path,
            "diamond": self.diamond_path,
            "samtools": self.samtools_path,
            "bedtools": self.bedtools_path,
            "psiblast": self.psiblast_path,
            "makeblastdb": self.makeblastdb_path,
        }.get(what, shutil.which(what))

    def __init__(
        self,
        *,
        mmseqs_path: Optional[str] = shutil.which("mmseqs"),
        diamond_path: Optional[str] = shutil.which("diamond"),
        blastn_path: Optional[str] = shutil.which("blastn"),
        blastp_path: Optional[str] = shutil.which("blastp"),
        blastx_path: Optional[str] = shutil.which("blastx"),
        tblastn_path: Optional[str] = shutil.which("tblastn"),
        tblastx_path: Optional[str] = shutil.which("tblastn"),
        psiblast_path: Optional[str] = shutil.which("psiblast"),
        samtools_path: Optional[str] = shutil.which("samtools"),
        bedtools_path: Optional[str] = shutil.which("bedtools"),
        makeblastdb_path: Optional[str] = shutil.which("makeblastdb"),
        ncbi_txdb_names_path: Optional[str] = None,
        ncbi_txdb_nodes_path: Optional[str] = None,
        chunk_size: int = CHUNK_SIZE_DEFAULT,
    ):
        self.mmseqs_path = mmseqs_path
        self.diamond_path = diamond_path
        self.blastn_path = blastn_path
        self.blastp_path = blastp_path
        self.blastx_path = blastx_path
        self.tblastn_path = tblastn_path
        self.tblastx_path = tblastx_path
        self.psiblast_path = psiblast_path
        self.samtools_path = samtools_path
        self.bedtools_path = bedtools_path
        self.makeblastdb_path = makeblastdb_path
        self.chunk_size = chunk_size
        self.ncbi_txdb_names_path = ncbi_txdb_names_path
        self.ncbi_txdb_nodes_path = ncbi_txdb_nodes_path


class AbstractBlastDB(abc.ABC):
    @abstractmethod
    def finalize(self):
        raise NotImplementedError

    @abstractmethod
    def search_exec(
        self,
        query_path: str,
        dest_path: str,
        bsc: BlastSearchConfig,
        tmpdir: Optional[str] = None,
    ):
        raise NotImplementedError

    @abstractmethod
    def gc(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def bc(self) -> BlastConfig:
        raise NotImplementedError
