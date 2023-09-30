#!/usr/bin/env python
import copy
import functools
import gc
import os
import queue
from collections import defaultdict

import msgpack
import pandas as pd
from tqdm import tqdm

from labw_utils.commonutils.lwio.safe_io import get_reader, get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import (
    Dict,
    Iterable,
    Tuple,
    Optional,
    Sequence,
    Mapping,
    Literal,
    Any,
)

_lh = get_logger(__name__)

IS_INHERITED = -1
NON_EXIST = -2
HOMO_SAPIENS = 9606
VERTEBRATA = 7742
EUKARYOTA = 2759
CELLULAR_ORGANISMS = 131567


NEW_TAXDB_FMT = {
    "nodes": [
        "tax_id",
        "parent_tax_id",
        "rank",
        "embl_code",
        "division_id",
        "inherited_div_flag",
        "genetic_code_id",
        "inherited_gc_flag",
        "mitochondrial_genetic_code_id",
        "inherited_mgc_flag",
        "genbank_hidden_flag",
        "hidden_subtree_root_flag",
        "comments",
        "plastid_genetic_code_id",
        "inherited_pgc_flag",
        "specified_species",
        "hydrogenosome_genetic_code_id",
        "inherited_hgc_flag",
    ],
    "names": ["tax_id", "name_txt", "unique_name", "name_class"],
    "division": ["division_id", "division_cde", "division_name", "comments"],
    "gencode": [
        "genetic code id",
        "abbreviation",
        "name",
        "cde",
        "starts",
    ],
    "delnodes": ["taxid"],
    "merged": ["old_tax_id", "new_tax_id"],
    "citations": [
        "cit_id",
        "cit_key",
        "medline_id",
        "pubmed_id",
        "url",
        "text",
        "taxid_list",
    ],
    "typeoftype": ["type_name", "synonyms", "nomenclature", "description"],
    "host": ["tax_id", "potential_hosts"],
    "typematerial": ["tax_id", "tax_name", "type", "identifier"],
    "rankedlineage": [
        "tax_id",
        "tax_name",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
        "superkingdom",
    ],
    "fullnamelineage": ["tax_id", "tax_name", "lineage"],
    "taxidlineage": ["tax_id", "lineage"],
    "excludedfromtype": ["tax_id", "tax_name", "property", "voucher_strain"],
}


class TaxonDB:
    _dfs: Dict

    def __init__(self, dfs: Dict, should_copy: bool = True):
        self._dfs = dfs if not should_copy else copy.deepcopy(dfs)

    def get_parent(self, txid: int) -> int:
        return self._node_txid_to_parent.get(txid, NON_EXIST)

    def get_child(self, txid: int) -> Sequence[int]:
        return self._node_txid_to_child.get(txid, [])

    def get_child_recursive(self, txid: int) -> Iterable[int]:
        cs = self.get_child(txid)
        for c in cs:
            yield c
            yield from self.get_child_recursive(c)

    def get_child_recursive_l(self, txid: int) -> Iterable[int]:
        q = queue.Queue()
        q.put(txid)
        while not q.empty():
            curr_txid = q.get()
            cs = self.get_child(curr_txid)
            for c in cs:
                q.put(c)
            yield curr_txid

    def get_division(self, txid: int) -> Optional[int]:
        return self._node_txid_to_division.get(txid)

    def get_genetic_code(self, txid: int) -> Optional[int]:
        return self._node_txid_to_genetic_code.get(txid)

    @property
    def _rankedlineage(
        self,
    ) -> Mapping[int, Tuple[str, str, str, str, str, str, str, str, str]]:
        """
        tax_id -> {tax_name, species, genus, family, order, class, phylum, kingdom, superkingdom}
        """
        return self._dfs["RL_TXID_TO_LINEAGE"]

    @property
    def _division(self) -> Mapping[int, Tuple[str, str]]:
        return self._dfs["DIVISION"]

    @property
    def _node_txid_to_parent(self) -> Mapping[int, int]:
        return self._dfs["NODE_TXID_TO_PARENT"]

    @property
    def _node_txid_to_child(self) -> Mapping[int, Sequence[int]]:
        return self._dfs["NODE_TXID_TO_CHILD"]

    @property
    def _node_txid_to_division(self) -> Mapping[int, int]:
        return self._dfs["NODE_TXID_TO_DIVISION"]

    @property
    def _node_txid_to_genetic_code(self) -> Mapping[int, int]:
        return self._dfs["NODE_TXID_TO_GENETIC_CODE"]

    @property
    def _node_division_to_txid(self) -> Mapping[int, Sequence[int]]:
        return self._dfs["NODE_DIVISION_TO_TXID"]

    @property
    def _node_txid_to_rank(self):
        return self._dfs["NODE_TXID_TO_RANK"]

    @property
    def _removed_merged(self) -> Mapping[int, int]:
        return self._dfs["REMOVED_MERGED"]

    @property
    def taxons(self) -> Sequence[int]:
        return self._dfs["NODE_TXID_TO_DIVISION"].keys()

    def limit_to(self, txids: Iterable[int]) -> Dict:
        txids = set(txids)
        new_dfs = {
            "REMOVED_MERGED": self._removed_merged,
            "NODE_TXID_TO_PARENT": {
                child_id: parent_id
                for child_id, parent_id in self._node_txid_to_parent.items()
                if parent_id in txids and child_id in txids
            },
            "RL_TXID_TO_LINEAGE": {k: v for k, v in self._rankedlineage.items() if k in txids},
            "NODE_TXID_TO_CHILD": {
                parent_id: list(child_id for child_id in child_ids if child_id in txids)
                for parent_id, child_ids in self._node_txid_to_child.items()
                if parent_id in txids
            },
            "NODE_TXID_TO_DIVISION": {k: v for k, v in self._node_txid_to_division.items() if k in txids},
            "NODE_TXID_TO_GENETIC_CODE": {k: v for k, v in self._node_txid_to_genetic_code.items() if k in txids},
            "NODE_TXID_TO_RANK": {k: v for k, v in self._node_txid_to_rank.items() if k in txids},
            "NODE_DIVISION_TO_TXID": defaultdict(lambda: []),
        }
        for k, v in new_dfs["NODE_TXID_TO_DIVISION"].items():
            new_dfs["NODE_DIVISION_TO_TXID"][v].append(k)
        new_dfs["NODE_DIVISION_TO_TXID"] = dict(new_dfs["NODE_DIVISION_TO_TXID"])
        new_dfs["DIVISION"] = {k: v for k, v in self._division.items() if k in new_dfs["NODE_DIVISION_TO_TXID"].keys()}
        return new_dfs

    def rebase(self, new_root_txid: int) -> "TaxonDB":
        _lh.info("REBASE: FETCHING CHILDREN...")
        affected_taxons = set(self.get_child_recursive(new_root_txid))
        affected_taxons.add(new_root_txid)
        affected_taxons.add(NON_EXIST)
        _lh.info("REBASE: LIMITING...")
        new_dfs = self.limit_to(affected_taxons)
        new_dfs["NODE_TXID_TO_PARENT"][new_root_txid] = new_root_txid
        new_instance = self.__class__(new_dfs, should_copy=False)
        _lh.info("REBASE: FINISHED")
        return new_instance

    def get_parent_traceback(self, txid: int) -> Iterable[int]:
        yield txid
        while True:
            parent_txid = self.get_parent(txid)
            if parent_txid is None or parent_txid == NON_EXIST or parent_txid == txid:
                return
            else:
                yield parent_txid
                txid = parent_txid

    # FIXME: have bugs!
    def find_root(self) -> Iterable[int]:
        for txid in tqdm(self.taxons, desc="Finding root..."):
            if self.get_parent(txid) == NON_EXIST or self.get_parent(txid) == txid:
                yield txid

    def trim_to(self, txids: Sequence[int]) -> "TaxonDB":
        """
        Keep only selected species and its parent.

        Remove irrelevant or child species.

        :param txids:
        :return:
        """
        _lh.info("TRIM TO: FETCHING PARENTS...")
        affected_taxons = set()
        txids = set(txids)
        for txid in tqdm(txids, "FETCHING PARENTS..."):
            if txid not in self.taxons:
                continue
            if txid not in affected_taxons:
                affected_taxons.update(self.get_parent_traceback(txid))
        affected_taxons.add(NON_EXIST)
        _lh.info("TRIM TO: LIMITING...")
        new_dfs = self.limit_to(affected_taxons)
        new_instance = self.__class__(new_dfs, should_copy=False)
        _lh.info("TRIM TO: FINISHED")
        return new_instance

    @classmethod
    def from_zsv(cls, dbpath: str):
        _lh.info("CONVERT DB -> MSGPACK")
        dfs = {}

        def load(dbname: str) -> pd.DataFrame:
            in_path = os.path.join(dbpath, f"{dbname}.dmp.zsv.gz")
            _lh.info("LOAD_DB '%s'", in_path)
            return pd.read_csv(
                in_path,
                sep="\0",
                names=NEW_TAXDB_FMT[dbname],
                engine="c",  # PyArrow does not support sep="\0"
            )

        def convert_removed_merged_table():
            _lh.info("CURRENT CONVERT: REMOVED_MERGED")
            removed_table = load("delnodes")
            dfs["REMOVED_MERGED"] = {}
            """txid -> [(name_txt, name_class)]"""
            for tup in tqdm(list(removed_table.itertuples()), "PARSING delnodes TABLE..."):
                dfs["REMOVED_MERGED"][tup.taxid] = NON_EXIST
            del removed_table
            gc.collect()

            merged_table = load("merged")
            for tup in tqdm(list(merged_table.itertuples()), "PARSING merged TABLE..."):
                dfs["REMOVED_MERGED"][tup.old_tax_id] = tup.new_tax_id
            del merged_table
            gc.collect()

        def convert_rankedlineage_table():
            _lh.info("CURRENT CONVERT: RL_TXID_TO_LINEAGE")
            rankedlineage_table = load("rankedlineage")
            dfs["RL_TXID_TO_LINEAGE"] = {NON_EXIST: tuple(["NON_EXIST"] * 9)}
            for tup in tqdm(
                list(rankedlineage_table.itertuples(index=False)),
                "PARSING rankedlineage TABLE...",
            ):
                dfs["RL_TXID_TO_LINEAGE"][tup.tax_id] = tuple(
                    map(lambda x: "UNNAMED" if x is pd.isna(x) else x, tup[1:])
                )
            del rankedlineage_table
            gc.collect()

        def convert_division_table():
            _lh.info("CURRENT CONVERT: DIVISION")
            division_table = load("division")
            dfs["DIVISION"] = {}
            for tup in tqdm(list(division_table.itertuples()), "PARSING division TABLE..."):
                dfs["DIVISION"][tup.division_id] = (tup.division_cde, tup.division_name)
            dfs["DIVISION"][NON_EXIST] = ("NAN", "NON_EXIST")
            del division_table
            gc.collect()

        def convert_node_table():
            _lh.info("CURRENT CONVERT: NODE")
            nodes_table = load("nodes")
            dfs["NODE_TXID_TO_PARENT"] = {NON_EXIST: NON_EXIST}
            dfs["NODE_TXID_TO_CHILD"] = defaultdict(lambda: [])
            dfs["NODE_TXID_TO_CHILD"][NON_EXIST] = []
            dfs["NODE_TXID_TO_DIVISION"] = {NON_EXIST: NON_EXIST}
            dfs["NODE_TXID_TO_RANK"] = {NON_EXIST: "NON_EXIST"}
            dfs["NODE_TXID_TO_GENETIC_CODE"] = {NON_EXIST: NON_EXIST}
            dfs["NODE_DIVISION_TO_TXID"] = defaultdict(lambda: [])
            dfs["NODE_DIVISION_TO_TXID"][NON_EXIST].append(NON_EXIST)

            for tup in tqdm(list(nodes_table.itertuples()), "PARSING nodes TABLE 1/3..."):
                dfs["NODE_TXID_TO_PARENT"][tup.tax_id] = tup.parent_tax_id
                dfs["NODE_TXID_TO_CHILD"][tup.parent_tax_id].append(tup.tax_id)
                dfs["NODE_TXID_TO_RANK"][tup.tax_id] = tup.rank

                if tup.division_id is None:
                    if tup.inherited_div_flag == 1:
                        this_division_id = IS_INHERITED
                    else:
                        this_division_id = None
                else:
                    this_division_id = tup.division_id
                dfs["NODE_TXID_TO_DIVISION"][tup.tax_id] = this_division_id

                if tup.genetic_code_id is None:
                    if tup.inherited_gc_flag == 1:
                        this_genetic_code_id = IS_INHERITED
                    else:
                        this_genetic_code_id = None
                else:
                    this_genetic_code_id = tup.genetic_code_id
                dfs["NODE_TXID_TO_GENETIC_CODE"][tup.tax_id] = this_genetic_code_id
            dfs["NODE_TXID_TO_CHILD"] = dict(dfs["NODE_TXID_TO_CHILD"])
            del nodes_table
            gc.collect()

            # Start resolving TXID_TO_DIVISION
            partial_instance = cls(dfs, should_copy=False)

            for txid in tqdm(partial_instance.taxons, "PARSING nodes TABLE 2/3..."):
                parent_txid = txid
                trace = [parent_txid]
                last_division = partial_instance.get_division(parent_txid)
                while last_division == IS_INHERITED:
                    parent_txid = partial_instance.get_parent(parent_txid)
                    trace.append(parent_txid)
                    if parent_txid is None:
                        parent_txid = 1
                        _lh.warning(
                            "%d failed in finding valid parent! TRACE: %s",
                            txid,
                            str(trace),
                        )
                    last_division = partial_instance.get_division(parent_txid)
                    if last_division != IS_INHERITED:
                        dfs["NODE_TXID_TO_DIVISION"][txid] = last_division
                dfs["NODE_DIVISION_TO_TXID"][last_division].append(txid)

            for txid in tqdm(partial_instance.taxons, "PARSING nodes TABLE 3/3..."):
                parent_txid = txid
                trace = [parent_txid]
                last_genetic_code = partial_instance.get_genetic_code(parent_txid)
                while last_genetic_code == IS_INHERITED:
                    parent_txid = partial_instance.get_parent(parent_txid)
                    trace.append(parent_txid)
                    if parent_txid is None:
                        _lh.warning(
                            "%d failed in finding valid parent! TRACE: %s",
                            txid,
                            str(trace),
                        )
                    last_genetic_code = partial_instance.get_genetic_code(parent_txid)
                    if last_genetic_code != IS_INHERITED:
                        dfs["NODE_TXID_TO_GENETIC_CODE"][txid] = last_genetic_code
            del partial_instance
            gc.collect()

        convert_rankedlineage_table()
        convert_division_table()
        convert_node_table()
        convert_removed_merged_table()

        # Finished
        return cls(dfs)

    def resolve_taxid(self, txid: int) -> int:
        """
        Resolve taxon ID that is merged or deletec.
        """
        resolved_id = self._removed_merged.get(txid, txid)
        if resolved_id not in self.taxons:
            resolved_id = NON_EXIST
        if resolved_id != txid:
            _lh.warning("TXID RESOLVE CHANGED %d -> %d", txid, resolved_id)
        return resolved_id

    @functools.lru_cache(maxsize=1024)
    def get_name(
        self,
        txid: int,
        lineage_level: Literal[
            "tax_name",
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "kingdom",
            "superkingdom",
        ] = "tax_name",
    ) -> Optional[str]:
        txid = self.resolve_taxid(txid)
        ranked_names = self._rankedlineage.get(txid)
        if ranked_names is None:
            _lh.warning(f"ERROR: {txid} NOT RECORDED IN RANKEDLINEAGE TABLE!")
            return None

        return ranked_names[
            [
                "tax_name",
                "species",
                "genus",
                "family",
                "order",
                "class",
                "phylum",
                "kingdom",
                "superkingdom",
            ].index(lineage_level)
        ]

    def to_msgpack(self, out_fn: str):
        _lh.info("DUMPING MSGPACK")
        with get_writer(out_fn, is_binary=True) as jwriter:
            msgpack.dump(self._dfs, jwriter)
        _lh.info("DUMPING MSGPACK: FINISHED")

    # def to_ete3(self, root_txid: int, root_txname: str) -> ete3.Tree:
    #     t = ete3.Tree()
    #     nodes = {}
    #     root = t.add_child(name="-".join((str(root_txid), root_txname)))
    #     nodes[root_txid] = root
    #     to_traverse = queue.Queue()
    #     pdb = tqdm(total=len(self), desc="Converting...")
    #     child_node_txids = list(self.get_child(txid=root_txid))
    #     for child_node_txid in child_node_txids:
    #         to_traverse.put(child_node_txid)
    #     while not to_traverse.empty():
    #         pdb.update(1)
    #         this_node_txid = to_traverse.get()
    #         # parent_node_id = self.get_parent(this_node_txid)
    #         # print(parent_node_id, this_node_txid)
    #         parent_node = nodes[self.get_parent(this_node_txid)]
    #         this_node_name = (
    #             self.get_name(this_node_txid)
    #             .replace("(", "_")
    #             .replace(")", "_")
    #             .replace("[", "_")
    #             .replace("]", "_")
    #             .replace("'", "_")
    #         )
    #         # print(this_node_name)
    #         nodes[this_node_txid] = parent_node.add_child(
    #             name="-".join((str(this_node_txid), this_node_name))
    #         )
    #         child_node_txids = self.get_child(txid=this_node_txid)
    #         for child_node_txid in child_node_txids:
    #             to_traverse.put(child_node_txid)
    #     pdb.close()
    #     return root

    @classmethod
    def from_msgpack(cls, in_fn: str):
        _lh.info("LOADING MSGPACK")
        with get_reader(in_fn, is_binary=True) as reader:
            new_dfs = msgpack.load(reader, strict_map_key=False)
        new_instance = cls(new_dfs, should_copy=False)
        _lh.info("LOADING MSGPACK: FINISHED")
        return new_instance

    def tree_view(self, root_txid: int) -> Mapping[Tuple[str, int], Any]:
        if root_txid == NON_EXIST:
            return {}
        retv = []
        for child_txid in self.get_child(root_txid):
            retv.append(self.tree_view(child_txid))

        return {(self.get_name(root_txid), root_txid): retv}

    def filter(self, txids: Iterable[int]) -> Iterable[int]:
        return filter(lambda x: x in self.taxons, txids)

    def __len__(self):
        return len(self.taxons)

    def __contains__(self, txid: int):
        return txid in self.taxons
