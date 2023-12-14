import os.path
import re
import time
from labw_utils.typing_importer import Optional

from Bio import GenBank
from Bio import Entrez
from Bio.GenBank.Record import Record as GBRecord, Feature

from blast_utils.blast6_io import parse_blast6, read_blast6
from blast_utils.blast_db import AbstractBlastDB, BlastConfig, BlastSearchConfig
from blast_utils.local_genbank_db import LocalGenBankDB
from blast_utils.ncbi_taxdb import NON_EXIST
from labw_utils.commonutils.lwio import get_reader
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import easyexec

_TXID_REGEX = re.compile(r'.*/db_xref="taxon:(\d+)".*')

Entrez.email = "theafamily@126.com"

_lh = get_logger(__name__)


class RemoteBlastDB(AbstractBlastDB):
    _lgdb: LocalGenBankDB
    _bc: BlastConfig

    def finalize(self):
        """Disabled for remote BLAST"""
        pass

    def __init__(self, local_genbank_db_path: str, bc: BlastConfig):
        self._lgdb = LocalGenBankDB(local_genbank_db_path)
        self._bc = bc

    def search_exec(
        self,
        query_path: str,
        dest_path: str,
        bsc: BlastSearchConfig,
        tmpdir: Optional[str] = None,
    ):
        _ = tmpdir
        del tmpdir
        blast_result_file = os.path.join(dest_path, "search_result.blast6.tsv")
        if bsc.search_executable not in {
            "blastp",
            "blastx",
            "blastn",
            "tblastn",
            "tblastx",
            # "psiblast", # Not supported.
        }:
            raise NotImplementedError("Operation not supported!")
        retv = 1
        n_retries = 0
        while not retv == 0 and n_retries < bsc.remote_blast_n_max_retries:
            n_retries += 1
            log_path = os.path.join(dest_path, f"{bsc.search_executable}_search.{n_retries}.log")
            _lh.info("Remote BLAST Attempt %d", n_retries)
            retv = easyexec(
                [
                    self._bc.which(bsc.search_executable),
                    "-query",
                    query_path,
                    "-db",
                    bsc.remote_blast_db_name,
                    "-out",
                    blast_result_file,
                    "-outfmt",
                    str(7),
                    "-evalue",
                    str(bsc.filter_e_value),
                    "-remote",
                    *bsc.cmdargs,
                ],
                log_path=log_path,
            )
            if retv != 0:
                time.sleep(10)
        result_blast6 = read_blast6(blast_result_file)
        self._lgdb.prepare(result_blast6["sseqid"])

        title_dict = {}
        txid_dict = {}
        seq_dict = {}
        sseqlen_dict = {}
        for tup in result_blast6.itertuples(index=False):
            gbacc = tup.sseqid
            with get_reader(self._lgdb.get(gbacc), is_binary=False) as reader:
                record: GBRecord = GenBank.read(reader)
                effective_txid = NON_EXIST
                for feature in record.features:
                    feature: Feature
                    if feature.key == "source":
                        for q in feature.qualifiers:
                            match_result = _TXID_REGEX.match(str(q))
                            if match_result is not None:
                                effective_txid = int(match_result.group(1))
                                break
            title_dict[gbacc] = record.definition
            txid_dict[gbacc] = effective_txid
            seq_dict[f"{gbacc}:{tup.sstart - 1}-{tup.send}"] = record.sequence[tup.sstart - 1 : tup.send]
            sseqlen_dict[gbacc] = len(record.sequence)

        parse_blast6(
            src_blast6_df=result_blast6,
            title_dict=title_dict,
            txid_dict=txid_dict,
            seq_dict=seq_dict,
            bsc=bsc,
        ).to_csv(
            os.path.join(dest_path, "search_result.parsed.tsv"),
            sep="\t",
            index=False,
        )

    def gc(self):
        """Disabled for remote BLAST"""
        pass

    @property
    def bc(self) -> BlastConfig:
        return self._bc
