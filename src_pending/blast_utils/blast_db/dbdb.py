from __future__ import annotations

import glob
import json
import multiprocessing
import os.path
import re
import socket
import subprocess
import threading
import time
from typing import Literal, Optional

from blast_utils import merge_table
from blast_utils.blast_db import BlastSearchConfig, BlastConfig, AbstractBlastDB
from blast_utils.blast_db._sbdb import SingleBlastDB
from labw_utils.bioutils.parser.fasta import FastaIterator
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.itertools_helper import window, head
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import parallel_map
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
from labw_utils.typing_importer import Iterable, Tuple, Dict

_lh = get_logger(__name__)

UNIPROT_TAXON_REGEX = re.compile(r".*OX=([0-9]+).*")
ENSEMBL_GENE_SYMBOL_REGEX = re.compile(r".*gene_symbol:(\S+).*")
ENSEMBL_GENE_ACCESSION_REGEX = re.compile(r".*gene:(\S+).*")
ENSEMBL_TRANSCRIPT_ACCESSION_REGEX = re.compile(r".*transcript:(\S+).*")


# def draw_taxa_tree_from_blast_hit(
#     result_tsv_path: str,
#     txdb_msgpack_path: str,
#     out_newick_path: str,
#     txdb_root: int = CELLULAR_ORGANISMS,
# ):
#     df = pd.read_csv(result_tsv_path, sep="\t")
#     queried_txid_set = set(df["txid"])
#     (
#         TaxonDB.from_msgpack(txdb_msgpack_path)
#         .rebase(txdb_root)
#         .trim_to(queried_txid_set)
#         .to_ete3(txdb_root, "root")
#         .write(format=8, outfile=out_newick_path)
#     )


class GetLastDBIDServer(threading.Thread):
    _fifo_path: str
    _split_db_path: str
    _is_ready: bool

    @property
    def is_ready(self) -> bool:
        return self._is_ready

    @property
    def fifo_path(self) -> str:
        return self._fifo_path

    @property
    def split_db_path(self) -> str:
        return self._split_db_path

    def generate(self):
        i = 0
        while os.path.exists(os.path.join(self.split_db_path, str(i))):
            i += 1
        os.makedirs(os.path.join(self.split_db_path, str(i)))
        _lh.info("Assigned DBID %d", i)
        return i

    def __init__(self, split_db_path: str, fifo_path: str):
        super().__init__()
        self._split_db_path = split_db_path
        self._fifo_path = fifo_path
        self._is_ready = False

    def run(self):
        sock = socket.socket(family=socket.AF_UNIX, type=socket.SOCK_STREAM)
        sock.bind(self._fifo_path)
        sock.listen(1)
        self._is_ready = True
        while True:
            conn, _ = sock.accept()
            msg_recv = json.loads(conn.recv(1024))
            if msg_recv["type"] == "TERM":
                return
            else:
                conn.send(bytes(json.dumps({"id": self.generate()}), encoding="UTF-8"))


def start_get_last_dbid_server(split_db_path: str, id_server_path: str):
    _lh.info("Starting GetLastDBIDServer")
    if os.path.exists(id_server_path):
        os.unlink(id_server_path)
    dbs = GetLastDBIDServer(split_db_path=split_db_path, fifo_path=id_server_path)
    dbs.start()
    while not dbs.is_ready:
        time.sleep(0.1)
    _lh.info("Started GetLastDBIDServer")
    return dbs


def stop_get_last_dbid_server(dbs: GetLastDBIDServer):
    _lh.info("Stopping GetLastDBIDServer")
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    sock.connect(dbs.fifo_path)
    sock.sendall(bytes(json.dumps({"type": "TERM"}), encoding="utf-8"))
    dbs.join()
    os.unlink(dbs.fifo_path)
    _lh.info("Stopped GetLastDBIDServer")


class DistributedBlastDB(AbstractBlastDB):
    _single_blast_dbs: Dict[int, SingleBlastDB]
    _split_db_path: str
    _bc: BlastConfig

    @property
    def split_db_path(self) -> str:
        return self.split_db_path

    @property
    def bc(self) -> BlastConfig:
        return self._bc

    def finalize(self, n_jobs: int = multiprocessing.cpu_count()):
        parallel_map(
            f=lambda a: a.finalize(),
            input_iterable=tqdm(self._single_blast_dbs.values(), desc="Finalizing..."),
            backend="threading",
            n_jobs=n_jobs,
        )

    def add_from_iterator(
        self,
        it: Iterable[Tuple[str, str, int, str]],
        id_server_path: str,
    ):
        """
        :param it: Iterable of [accession, title, txid, seq]
        """
        _lh.info("Adding from iterator with chunk size %d", self._bc.chunk_size)
        pb = tqdm(total=None, desc="Splitting...")
        for ls in window(it, size=self._bc.chunk_size, last_action="ignore"):
            sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            sock.connect(id_server_path)
            sock.sendall(bytes(json.dumps({"type": "REQ"}), encoding="utf-8"))
            dbid = json.loads(sock.recv(1024))["id"]
            sock.close()

            self._single_blast_dbs[dbid] = SingleBlastDB(
                tax_split_path=self._split_db_path,
                dbid=dbid,
                db_type=self._db_type,
                bc=self._bc,
            )
            for l in ls:
                accession, title, txid, seq = l
                self._single_blast_dbs[dbid].add(
                    accession=accession,
                    title=title,
                    txid=txid,
                    seq=seq,
                )
                pb.update(1)

    def add_from_ensembl_reference(
        self,
        ensembl_fa_path: str,
        txid: int,
        taxon_name: str,
        id_server_path: str,
        first_n_records_only: int = 0,
    ):
        def generate_sym_acc(_seqid: str):
            transcript_acc_match = ENSEMBL_TRANSCRIPT_ACCESSION_REGEX.match(_seqid)
            gene_acc_match = ENSEMBL_GENE_ACCESSION_REGEX.match(_seqid)
            gene_sym_match = ENSEMBL_GENE_SYMBOL_REGEX.match(_seqid)
            real_accession = _seqid.split(" ")[0]
            rets = "-".join(
                [
                    gene_sym_match.group(1) if gene_sym_match is not None else "UK_GENE_SYM",
                    gene_acc_match.group(1) if gene_acc_match is not None else "UK_GENE_ACC",
                    transcript_acc_match.group(1) if transcript_acc_match is not None else "UK_TRANS_ACC",
                    real_accession,
                ]
            )
            return rets

        def get_iterable() -> Iterable[Tuple[str, str, int, str]]:
            with FastaIterator(ensembl_fa_path, show_tqdm=False) as fai:
                for record in fai:
                    accession = record.seq_id.split(" ")[0]
                    title = f"{taxon_name}-{generate_sym_acc(record._seq_id)}"
                    yield accession, title, txid, record.sequence

        _lh.info("Adding ENSEMBL reference from %s", ensembl_fa_path)
        it = get_iterable()
        if first_n_records_only != 0:
            it = head(it, n=first_n_records_only)
        self.add_from_iterator(it, id_server_path)
        _lh.info("Adding ENSEMBL reference from %s FIN", ensembl_fa_path)

    def _add_from_prebuilt_blastdb(
        self,
        blastdb_path: str,
        first_n_records_only: int,
        id_server_path: str,
    ):
        def lex_record(s: bytes) -> Tuple[str, str, int, str]:
            try:
                s = s.decode(encoding="utf-8")
            except UnicodeDecodeError:
                s_err = s.decode(encoding="utf-8", errors="ignore")
                _lh.error(f"Error in decoding bytes: '{s_err}'. This record is ignored.")
            x = s.strip().split("\t")
            return x[0], x[1], int(x[2]), x[3]

        _lh.info("Adding BLAST-DB reference from %s", blastdb_path)
        p = subprocess.Popen(
            [
                "blastdbcmd",
                "-db",
                blastdb_path,
                "-outfmt",
                "\t".join(
                    (
                        "%a",  # Accession
                        "%t",  # Title
                        "%T",  # TaxID
                        "%s",  # Sequence
                    )
                ),
                "-entry",
                "all",
            ],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if first_n_records_only != 0:
            it = map(lex_record, head(p.stdout, n=first_n_records_only))
        else:
            it = map(lex_record, p.stdout)
        self.add_from_iterator(it, id_server_path)
        if first_n_records_only != 0:
            p.kill()
            p.wait()
        else:
            if p.wait() != 0:
                raise RuntimeError(
                    "Failed read to BLAST database! Details:\n"
                    + str(p.stderr.read(), encoding="utf-8", errors="ignore")
                )
        _lh.info("Adding BLAST-DB reference from %s FIN", blastdb_path)

    def add_from_prebuilt_blastdb(
        self,
        blastdb_path: str,
        id_server_path: str,
        first_n_records_only: int = 0,
    ):
        if first_n_records_only != 0:
            self._add_from_prebuilt_blastdb(
                blastdb_path=blastdb_path,
                first_n_records_only=first_n_records_only,
                id_server_path=id_server_path,
            )
        else:
            p = subprocess.Popen(
                [
                    "blastdbcmd",
                    "-db",
                    blastdb_path,
                    "-metadata",
                ],
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if p.wait() != 0:
                raise RuntimeError(
                    "Failed read to BLAST database! Details:\n"
                    + str(p.stderr.read(), encoding="utf-8", errors="ignore")
                )
            else:
                try:
                    db_metadata = json.load(p.stdout)
                    n_volumes = db_metadata["number-of-volumes"]
                except Exception:
                    raise RuntimeError("Failed read to BLAST database! Details: metadata decode failed")
                if n_volumes == 1:
                    _lh.info("Single-volume DB")
                    self._add_from_prebuilt_blastdb(
                        blastdb_path=blastdb_path,
                        first_n_records_only=0,
                        id_server_path=id_server_path,
                    )
                else:
                    _lh.info("Multi-volume (%d) DB", n_volumes)
                    n_digits = len(str(n_volumes))
                    parallel_map(
                        lambda i: self._add_from_prebuilt_blastdb(
                            blastdb_path=".".join((blastdb_path, str(i).rjust(n_digits, "0"))),
                            first_n_records_only=0,
                            id_server_path=id_server_path,
                        ),
                        range(n_volumes),
                        backend="loky",
                    )

    @classmethod
    def new(cls, split_db_path: str, db_type: Literal["prot", "nucl"], bc: BlastConfig):
        _lh.info("Removing previous database...")
        rm_rf(split_db_path)
        _lh.info("Creating new database...")
        return cls(split_db_path=split_db_path, db_type=db_type, bc=bc)

    def __init__(
        self,
        split_db_path: str,
        db_type: Literal["prot", "nucl"],
        bc: BlastConfig,
    ):
        self._split_db_path = os.path.abspath(split_db_path)
        os.makedirs(self._split_db_path, exist_ok=True)
        self._single_blast_dbs = {}
        self._db_type = db_type
        for p in glob.glob(os.path.join(self._split_db_path, "*")):
            dbid = int(os.path.basename(p))
            self._single_blast_dbs[dbid] = SingleBlastDB(
                tax_split_path=self._split_db_path,
                dbid=dbid,
                db_type=self._db_type,
                bc=bc,
            )
        self._bc = bc

    def search_exec(
        self,
        query_path: str,
        dest_path: str,
        bsc: BlastSearchConfig,
        tmpdir: Optional[str] = None,
    ):
        if tmpdir is None:
            tmpdir = os.path.join(dest_path, "tmp")
        os.makedirs(tmpdir, exist_ok=True)

        _lh.info("Start searching...")
        bsc_non_parallel = bsc.duplicate()
        bsc_non_parallel.n_jobs = 1
        parallel_map(
            f=lambda a: a.search_exec(
                bsc=bsc_non_parallel,
                query_path=query_path,
                dest_path=os.path.join(tmpdir, str(a.dbid)),
            ),
            input_iterable=tqdm(self._single_blast_dbs.values(), desc="Searching..."),
            backend="threading",
            n_jobs=bsc.n_jobs,
        )
        merge_table(
            os.path.join(tmpdir, "*", "search_result.parsed.tsv"),
            os.path.join(dest_path, "search_result.parsed.tsv"),
        )

    def gc(self):
        for single_blast_db in self._single_blast_dbs.values():
            single_blast_db.gc()
