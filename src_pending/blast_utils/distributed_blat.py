import enum
import glob
import multiprocessing
import os.path
import shutil
import uuid

import msgpack
import pandas as pd

from blast_utils import AbstractConfig, merge_table
from blast_utils.blast6_io import read_blast6
from blast_utils.ncbi_taxdb import NON_EXIST
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer, get_reader
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import parallel_map, easyexec
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
from labw_utils.typing_importer import List, Tuple, Dict, Mapping, Sequence, Optional

_lh = get_logger(__name__)

CHUNK_SIZE_DEFAULT = 1 << 20


def parse_blat_blast6(
    t: pd.DataFrame,
    contig_dict: Mapping[str, str],
    txid_dict: Mapping[str, int],
):
    t["qstart"] = t["qstart"] - 1
    t["sstart"] = t["sstart"] - 1
    t["contig_id"] = list(
        map(
            lambda x: contig_dict.get(x, "chrUn"),
            t["sseqid"],
        )
    )
    t["txid"] = list(map(lambda x: txid_dict.get(x, NON_EXIST), t["sseqid"]))
    return t


class BlatConfig(AbstractConfig):
    chunk_size: int
    blat_path: Optional[str]
    fa_to_2bit_path: Optional[str]

    def which(self, what: str) -> Optional[str]:
        return {
            "blat": self.blat_path,
            "faToTwoBit": self.fa_to_2bit_path,
        }.get(what, shutil.which(what))

    def __init__(
        self,
        chunk_size: int = CHUNK_SIZE_DEFAULT,
        blat_path: Optional[str] = shutil.which("blat"),
        fa_to_2bit_path: Optional[str] = shutil.which("faToTwoBit"),
    ):
        self.chunk_size = chunk_size
        self.blat_path = blat_path
        self.fa_to_2bit_path = fa_to_2bit_path


class ReferenceState(enum.Enum):
    pending = 1
    as_fasta = 2
    as_2bit = 3


class _SingleBlatDB:
    _bc: BlatConfig
    _path: str
    _dbid: int

    _pending: Dict[str, Tuple[str, str, int]]
    _reference_filename: str
    _reference_state: ReferenceState

    @property
    def dbid(self) -> int:
        return self._dbid

    def __init__(self, tax_split_path: str, dbid: int, bc: BlatConfig):
        self._dbid = dbid
        self._bc = bc
        self._path = os.path.join(tax_split_path, str(self._dbid))
        os.makedirs(self._path, exist_ok=True)

        self._reference_filename = ""
        self._reference_state = ReferenceState.pending

    def add(self, fa_filename: str, contig_name: str, txid: int):
        self._pending[str(uuid.uuid4())] = (fa_filename, contig_name, txid)

    def finalize(self):
        # Construct FASTA
        final_fa_filename = os.path.join(self._path, "seqs.fa")
        contig_dict = {}
        txid_dict = {}

        with FastaWriter(final_fa_filename) as faw:
            for contig_uuid, (
                fa_filename,
                contig_name,
                txid,
            ) in self._pending.items():
                with FastaViewFactory(fa_filename, read_into_memory=False, show_tqdm=False) as fav:
                    faw.write(
                        FastaRecord(
                            seq_id=contig_uuid,
                            sequence=fav.sequence(contig_name),
                        )
                    )
                contig_dict[contig_uuid] = contig_name
                txid_dict[contig_uuid] = txid
        self._reference_state = ReferenceState.as_fasta
        self._pending.clear()

        with get_writer(os.path.join(self._path, "contig.msgpack"), is_binary=True) as w:
            msgpack.dump(contig_dict, w)

        with get_writer(os.path.join(self._path, "txid.msgpack"), is_binary=True) as w:
            msgpack.dump(txid_dict, w)

        # Convert to 2bit
        if self._bc.fa_to_2bit_path is not None:
            twobit_filename = os.path.join(self._path, "seqs.2bit")
            easyexec(
                [self._bc.blat_path, "-long", final_fa_filename, twobit_filename],
                log_path=os.path.join(self._path, "faToTwoBit.log"),
            )
            self._reference_state = ReferenceState.as_2bit
        else:
            _lh.warning("faToTwoBit not found; build FASTA reference")
            self._reference_filename = final_fa_filename
        with get_writer(os.path.join(self._path, "metadata.msgpack"), is_binary=True) as w:
            msgpack.dump(self._pending, w)

    def search(
        self,
        query_path: str,
        dest_path: str,
        cmdargs: List[str],
    ):
        self._bc.require(["blat"])
        blat_result_file = os.path.join(dest_path, "search_result.blast6")
        easyexec(
            [
                self._bc.blat_path,
                self._reference_filename,
                query_path,
                blat_result_file,
                "-out=blast9",
                *cmdargs,
            ],
            log_path=os.path.join(dest_path, "blat_search.log"),
        )
        with get_reader(os.path.join(self._path, "contig.msgpack"), is_binary=True) as r:
            contig_dict = msgpack.load(r)
        with get_reader(os.path.join(self._path, "txid.msgpack"), is_binary=True) as r:
            txid_dict = msgpack.load(r)

        t = parse_blat_blast6(read_blast6(blat_result_file), contig_dict, txid_dict)
        t.to_csv(
            os.path.join(dest_path, "search_result.parsed.tsv"),
            sep="\t",
            index=False,
        )


class DistributedBlatDB:
    _single_blast_dbs: Dict[int, _SingleBlatDB]
    _split_db_path: str
    _bc: BlatConfig

    def _get_last_db_id(self):
        i = 0
        while os.path.exists(os.path.join(self._split_db_path, str(i))):
            i += 1
        return i

    def add(self, fasta_path: str, txid: int):
        self.add_multi([fasta_path], [txid])

    def add_multi(self, fasta_paths: Sequence[str], txid: Sequence[int]):
        assert len(fasta_paths) == len(txid)
        # TODO

    @classmethod
    def new(cls, split_db_path: str, bc: BlatConfig):
        rm_rf(split_db_path)
        return cls(split_db_path=split_db_path, bc=bc)

    def __init__(
        self,
        split_db_path: str,
        bc: BlatConfig,
    ):
        self._split_db_path = os.path.abspath(split_db_path)
        os.makedirs(self._split_db_path, exist_ok=True)
        self._single_blast_dbs = {}
        for p in glob.glob(os.path.join(self._split_db_path, "*")):
            dbid = int(os.path.basename(p))
            self._single_blast_dbs[dbid] = _SingleBlatDB(
                tax_split_path=self._split_db_path,
                dbid=dbid,
                bc=bc,
            )
        self._bc = bc

    def finalize(self, n_jobs: int = multiprocessing.cpu_count()):
        parallel_map(
            f=lambda a: a.finalize(),
            input_iterable=tqdm(self._single_blast_dbs.values(), desc="Finalizing..."),
            backend="threading",
            n_jobs=n_jobs,
        )

    def search(
        self,
        query_path: str,
        dest_path: str,
        cmdargs: List[str],
        n_jobs: int = multiprocessing.cpu_count(),
    ):
        tmpdir = os.path.join(dest_path, "tmp")
        os.makedirs(tmpdir, exist_ok=True)

        _lh.info("Start searching...")
        parallel_map(
            f=lambda a: a.search(
                query_path=query_path,
                dest_path=os.path.join(tmpdir, str(a.dbid)),
                cmdargs=cmdargs,
            ),
            input_iterable=tqdm(self._single_blast_dbs.values(), desc="Searching..."),
            backend="threading",
            n_jobs=n_jobs,
        )
        merge_table(
            os.path.join(tmpdir, "*", "search_result.parsed.tsv"),
            os.path.join(dest_path, "search_result.tsv"),
        )
