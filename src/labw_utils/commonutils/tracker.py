import itertools
import os.path
from hashlib import blake2b
from json import JSONDecodeError

from labw_utils import Dict
from labw_utils.commonutils.lwio import get_reader
from labw_utils.commonutils.serializer.json import AbstractJSONSerializable
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Sequence, Mapping, Any, List, Tuple, Optional

CHUNK_SIZE = 1 << 14

_lh = get_logger(__name__)


class FTracker(AbstractJSONSerializable):
    @staticmethod
    def _dump_versions() -> Optional[Mapping[str, Any]]:
        return {}

    @staticmethod
    def _dump_metadata() -> Optional[Mapping[str, Any]]:
        pass

    @staticmethod
    def _validate_versions(versions: Mapping[str, Any]) -> None:
        return

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "file_hash_map": self.file_hash_map,
            "source_dest_map": self.source_dest_map,
            "fgrp": self.fgrp,
        }

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

    _title = "FTracker"
    file_hash_map: Dict[str, str]
    source_dest_map: List[Tuple[str, str]]

    fgrp: Dict[str, Tuple[str]]

    def __init__(
        self,
        file_hash_map: Mapping[str, str],
        source_dest_map: Sequence[Tuple[str, str]],
        fgrp: Mapping[str, Sequence[str]],
        **kwargs,
    ):
        _ = kwargs
        del kwargs
        self.source_dest_map = list(tuple(x) for x in source_dest_map)
        self.fgrp = {k: tuple(v) for k, v in fgrp.items()}
        self.file_hash_map = dict(file_hash_map)

    @property
    def rev_fgrp(self) -> Mapping[Tuple[str], str]:
        return {v: k for k, v in self.fgrp.items()}

    @staticmethod
    def get_hash(src_fpath: str) -> str:
        this_hash = blake2b()
        with get_reader(src_fpath, is_binary=True) as r:
            while True:
                c = r.read(CHUNK_SIZE)
                if not c:
                    break
                this_hash.update(c)
        return this_hash.hexdigest()

    def is_up_to_date(self, source: Sequence[str], dest: Sequence[str]) -> Tuple[bool, str]:
        dest_fgid = self.rev_fgrp.get(tuple(sorted(dest)))
        if dest_fgid is None:
            return False, "Dest file group not recognized"
        source_fgid = self.rev_fgrp.get(tuple(sorted(source)))
        if source_fgid is None:
            return False, "Source files not recognized"
        if (source_fgid, dest_fgid) not in self.source_dest_map:
            return False, "Relation not recognized"

        for fn in itertools.chain(source, dest):
            if not os.path.exists(fn):
                return False, f"File {fn} not exist"
            if fn not in self.file_hash_map:
                return False, f"File {fn} not recognized"
            if FTracker.get_hash(fn) != self.file_hash_map[fn]:
                return False, f"File {fn} hash mismatch"
        return True, ""

    def add(self, source: Sequence[str], dest: Sequence[str]):
        source_fg = tuple(sorted(source))

        if self.rev_fgrp.get(source_fg) is None:
            self.fgrp[str(len(self.fgrp))] = source_fg
        source_fgid = self.rev_fgrp.get(source_fg)

        dest_fg = tuple(sorted(dest))
        if self.rev_fgrp.get(dest_fg) is None:
            self.fgrp[str(len(self.fgrp))] = dest_fg
        dest_fgid = self.rev_fgrp.get(dest_fg)

        if not (source_fgid, dest_fgid) in self.source_dest_map:
            self.source_dest_map.append((source_fgid, dest_fgid))
        for fn in itertools.chain(source, dest):
            self.file_hash_map[fn] = FTracker.get_hash(fn)


class FTrackerDaemon:
    ft: FTracker
    fn: str

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __init__(self, fn: str):
        self.fn = fn
        try:
            self.ft = FTracker.load(fn)
            _lh.info("FTrack recovered from %s", fn)
        except (FileNotFoundError, JSONDecodeError, TypeError) as e:
            _lh.info("FTrack recovered from %s FAILED due to %s", fn, e)
            self.ft = FTracker({}, [], {})
        self.flush()

    def flush(self):
        self.ft.save(self.fn)

    def close(self):
        self.flush()
