"""
Transposon database that manages conserved transposon sequences.

.. versionadded:: 3.2.0
"""

import json
import random
import re

from labw_utils import UnmetDependenciesError
from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio import get_writer, get_reader
from labw_utils.commonutils.lwio.safe_io import get_reader
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Mapping, Tuple, Optional, List, Collection

_lh = get_logger(__name__)


class DfamEmblRecord:
    embl_id: Optional[str]
    rmsk_name: Optional[str]
    rmsk_type: Optional[str]
    rmsk_subtype: Optional[str]
    rmsk_species: Optional[str]

    def __init__(self):
        self.embl_id = ""
        self.rmsk_name = ""
        self.rmsk_type = ""
        self.rmsk_subtype = ""
        self.rmsk_species = ""


class TransposonDatabase:
    _accession_sequence_map: Mapping[str, str]
    _accession_hmm_map: Mapping[str, str]
    _accessions: List[str]
    _hmm_epool: List[Tuple[str, str]]

    @staticmethod
    def convert_dfam_embl(
        *,
        src_dfam_embl_file_path: str,
        required_taxa: Collection[str],
        dst_index_file_path: str,
        dst_consensus_fa_path: Optional[str],
        with_tqdm: bool = True,
    ):
        _lh.info("Enumerating accessions...")
        accession_info = {}
        current_record = DfamEmblRecord()
        sequence = ""
        if with_tqdm:
            f = get_reader(src_dfam_embl_file_path)
        else:
            f = get_tqdm_line_reader(src_dfam_embl_file_path)

        rmsk_annotations = False
        for line in f:
            line: str
            line_info = line[2:].strip()
            line_identifier = line[:2]

            if line_identifier == "ID":
                current_record.embl_id = line.split(";")[0].split(" ")[0]
            elif line_identifier == "CC":
                if "RepeatMasker Annotations:" in line_info:
                    rmsk_annotations = True
                elif rmsk_annotations:
                    parts = line_info.split(":")
                    if len(parts) > 1:
                        value = parts[1].strip()
                        key = parts[0].strip().lower()
                        if key == "type":
                            current_record.rmsk_type = value
                        elif key == "subtype":
                            current_record.rmsk_subtype = value
                        elif key == "species":
                            current_record.rmsk_species = value
                        else:
                            rmsk_annotations = False
                else:
                    rmsk_annotations = False
            elif line_identifier == "NM":
                current_record.rmsk_name = line_info
            elif line_identifier == "  ":
                sequence += re.sub(r"\d", "", line.replace(" ", ""))
            elif line_identifier == "//":
                if (
                    current_record.embl_id != ""
                    and current_record.rmsk_species != ""
                    and sequence != ""
                    and current_record.rmsk_species in required_taxa
                ):
                    rmsk_name = (
                        current_record.rmsk_name if current_record.rmsk_name != "" else f"NAME-{current_record.embl_id}"
                    )
                    rmsk_type = current_record.rmsk_type if current_record.rmsk_type != "" else "UNKNOWN"
                    rmsk_subtype = (
                        current_record.rmsk_subtype if current_record.rmsk_subtype != "" else f"{rmsk_type}-UNKNOWN"
                    )
                    accession_info[rmsk_name] = {
                        "CONSENSUS": sequence,
                        "HMM": "",
                        "NAME": rmsk_name,
                        "TYPE": rmsk_type,
                        "SUBTYPE": rmsk_subtype,
                    }

                current_record = DfamEmblRecord()
                sequence = ""
        _lh.info(
            "Finished with %d accesions, writing...",
            len(accession_info),
        )
        with get_writer(dst_index_file_path, is_binary=False) as w:
            json.dump(accession_info, w)
        if dst_consensus_fa_path is not None:
            with FastaWriter(dst_consensus_fa_path) as faw:
                for v in accession_info.values():
                    faw.write(FastaRecord(f"{v['NAME']}#{v['TYPE']}/{v['SUBTYPE']}", v["CONSENSUS"]))
        _lh.info("Finished")

    @staticmethod
    def convert_dfam_hdf5(
        *,
        src_dfam_hdf5_file_path: str,
        required_txid: int,
        dst_index_file_path: str,
        dst_consensus_fa_path: Optional[str],
        dst_hmm_path: Optional[str],
        with_tqdm: bool = True,
        fetch_parent: bool = True,
    ) -> None:
        try:
            import h5py
        except ImportError as e:
            raise UnmetDependenciesError("h5py") from e
        _lh.info("Enumerating taxons...")
        txids = [required_txid]
        accession_info = {}
        with h5py.File(src_dfam_hdf5_file_path) as ds:
            while fetch_parent:
                this_taxon_node = ds["Taxonomy"]["Nodes"][str(txids[-1])]
                if "Parent" in this_taxon_node:
                    txids.append(this_taxon_node["Parent"][0])
                else:
                    break

            _lh.info("Enumerating taxons finished with %d taxons in list", len(txids))
            _lh.info("Enumerating accessions...")
            accessions = set()
            for txid in txids:
                this_taxon_node = ds["Taxonomy"]["Nodes"][str(txid)]
                if "Families" in this_taxon_node:
                    accessions.update(this_taxon_node["Families"])
            if with_tqdm:
                it = tqdm(accessions)
            else:
                it = accessions
            for accession in it:
                try:
                    d = ds["Families"][accession[0:2]][accession[2:4]][accession[4:6]][accession]
                except KeyError:
                    try:
                        d = ds["Families"]["ByName"][accession]
                    except KeyError:
                        _lh.warning(
                            "Accession '%s' neither resolved as name nor accession, skipped",
                            accession,
                        )
                        continue
                try:
                    rmsk_name = d.attrs.get("name", f"NAME-{accession}")
                    rmsk_type = d.attrs.get("repeat_type", "UNKNOWN")
                    rmsk_subtype = d.attrs.get("repeat_subtype", f"{rmsk_type}-UNKNOWN")
                    accession_info[rmsk_name] = {
                        "CONSENSUS": d.attrs["consensus"],
                        "HMM": d.attrs.get("model", ""),
                        "NAME": rmsk_name,
                        "TYPE": rmsk_type,
                        "SUBTYPE": rmsk_subtype,
                    }
                except KeyError:
                    _lh.warning(
                        "Name '%s' do not have consensus/model, skipped",
                        accession,
                    )
                    continue
            _lh.info(
                "Finished with %d accesions, writing...",
                len(accession_info),
            )
            with get_writer(dst_index_file_path, is_binary=False) as w:
                json.dump(accession_info, w)
            if dst_consensus_fa_path is not None:
                with FastaWriter(dst_consensus_fa_path) as faw:
                    for v in accession_info.values():
                        faw.write(FastaRecord(f"{v['NAME']}#{v['TYPE']}/{v['SUBTYPE']}", v["CONSENSUS"]))
            if dst_hmm_path is not None:
                with get_writer(dst_hmm_path, is_binary=False) as hmmw:
                    for v in accession_info.values():
                        hmmw.write(v["HMM"])
            _lh.info("Finished")

    @classmethod
    def load(cls, index_path: str):
        with get_reader(index_path, is_binary=False) as r:
            accession_info = json.load(r)

        return cls(
            {k: v["CONSENSUS"] for k, v in accession_info.items()},
            {k: v["HMM"] for k, v in accession_info.items()},
        )

    def dump_json(self, dst_json_file_path: str):
        with get_writer(dst_json_file_path, is_binary=False) as w:
            json.dump(self._accession_sequence_map, w, indent=4)

    def __init__(
        self,
        accession_sequence_map: Mapping[str, str],
        accession_hmm_map: Mapping[str, str],
        **kwargs,
    ) -> None:
        _ = kwargs
        del kwargs
        self._accession_sequence_map = accession_sequence_map
        self._accession_hmm_map = accession_hmm_map
        self._accessions = list(self._accession_sequence_map.keys())
        self._hmm_epool = []

    def draw(self) -> Tuple[str, str]:
        """
        Randomly pick one transposon from the database.

        :returns: Transposon accession, start position and its sequence.
        """
        rdg = random.SystemRandom()
        selected_accn = rdg.choice(self._accessions)
        selected_seq = self._accession_sequence_map[selected_accn]
        return selected_accn, selected_seq

    def seq(self, transposon_accession: str) -> str:
        return self._accession_sequence_map[transposon_accession]
