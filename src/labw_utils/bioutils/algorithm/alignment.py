"""
``labw_utils.bioutils.algorithm.alignment`` -- Python-Implemented Alignment Algorithms

Some Python-Implemented alignment functions, including Smith-Waterman,
with other utilities like edit distance

.. versionadded:: 1.0.2
"""

import functools
import glob
import math
import os.path

from labw_utils import UnmetDependenciesError
from labw_utils.typing_importer import Optional, List, Tuple, Dict, Iterable

try:
    import numpy as np
    import numpy.typing as npt
except ImportError:
    raise UnmetDependenciesError("numpy")

from labw_utils.commonutils.lwio.safe_io import get_reader


class SubstMatrix:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    _real_mtx: Dict[str, Dict[str, int]]

    def __init__(self, real_mtx: Dict[str, Dict[str, int]]):
        self._real_mtx = real_mtx

    def get(self, seq1: str, seq2: str) -> int:
        return self._real_mtx[seq1][seq2]

    @classmethod
    def parse(cls, mat_file_path: str):
        real_mtx = {}
        header = []
        with get_reader(mat_file_path) as reader:
            for line in reader:
                line = line.strip()
                if line.startswith("#") or line == "":
                    continue
                else:
                    if not header:
                        header = line.split(" ")
                        while "" in header:
                            header.remove("")
                    else:
                        contents = line.split(" ")
                        while "" in contents:
                            contents.remove("")
                        h1 = contents.pop(0)
                        real_mtx[h1] = {h2: v for h2, v in zip(header, map(int, contents))}
        return cls(real_mtx)


_FILE_DIR = os.path.dirname(os.path.abspath(__file__))

_SUBST_MTX_PATHS = {
    os.path.basename(fn).split(".")[0]: fn for fn in glob.glob(os.path.join(_FILE_DIR, "scoring_mtx", "*.mat"))
}


def get_subst_mtx_names() -> Iterable[str]:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    return iter(_SUBST_MTX_PATHS.keys())


@functools.lru_cache()
def get_subst_mtx(mtx_name: str) -> SubstMatrix:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    return SubstMatrix.parse(_SUBST_MTX_PATHS[mtx_name])


class SmithWatermanAligner:
    """
    A General-Purposed Smith-Waterman Aligner.

    >>> SmithWatermanAligner('AAA', 'AAA').get_backtrack()[0]
    '>aln:seq1:qual:seq2:15\\nAAA\\n===\\nAAA'
    >>> SmithWatermanAligner('AAA', 'ATA').get_backtrack()[0]
    '>aln:seq1:qual:seq2:6\\nA--AA\\n=IIDD\\nATA--'
    >>> SmithWatermanAligner('AAA', 'AA').get_backtrack()[0]
    '>aln:seq1:qual:seq2:10\\nAAA\\n==D\\nAA-'
    >>> SmithWatermanAligner('AA', 'TTTA').get_backtrack()[0]
    '>aln:seq1:qual:seq2:5\\n---AA\\nIII=D\\nTTTA-'
    >>> SmithWatermanAligner('JHBGTYBYTJAA', 'YJGYJTVJYAA').get_backtrack()[0]
    '>aln:seq1:qual:seq2:10\\n-JHBGTYBY-T-J-AA\\nI=DD=D=DDI=I=I==\\nYJ--G-Y--JTVJYAA'
    >>> SmithWatermanAligner('TATATATGCGGGTAATTTAGGGCGGATCATGA', 'ATGCGGC').get_backtrack()[0]
    '>aln:seq1:qual:seq2:30\\nTATATATGCGGGTAATTTAGGGCGGATCATGA\\nD==DDDD====MDDDDDDDDDDDDDDDDDDDD\\n-AT----GCGGC--------------------'

    .. versionadded:: 1.0.2
    """

    seq1: str
    """Sequence 1"""

    seq2: str
    """Sequence 2"""

    match_score: int
    """Score for match"""

    mismatch_score: int
    """Score for mismatch, should be negative"""

    indel_score: int
    """Score for insertions and deletions"""

    is_global: bool
    """Whether the alignment should be global"""

    _sw_matrix: Optional[npt.NDArray[int]]
    """The Smith-Waterman Alignment Matrix"""

    _score: Optional[int]
    """Final Smith-Waterman Alignment Score"""

    _backtrack: Optional[List[str]]
    """Final Smith-Waterman backtrack"""

    _subst_mtx: Optional[SubstMatrix]
    """Substitution Matrix"""

    def __init__(
        self,
        seq1: str,
        seq2: str,
        match_score: int = 5,
        mismatch_score: int = -4,
        indel_score: int = -4,
        subst_mtx: Optional[SubstMatrix] = None,
        is_global: bool = True,
    ):
        self.seq1 = seq1
        self.seq2 = seq2
        self.mismatch_score = mismatch_score
        self.match_score = match_score
        self.indel_score = indel_score
        self.is_global = is_global
        self._sw_matrix = None
        self._score = None
        self._backtrack = None
        self._subst_mtx = subst_mtx

    def alignment(self) -> int:
        """
        Perform alignment and get the score
        """
        if self._sw_matrix is None:
            self._build_smith_waterman_matrix()
        self._score = np.max(self._sw_matrix[1:, 1:])
        return self._score

    @property
    def score(self) -> int:
        """
        Get Smith-Waterman Score.
        If alignment is not performed, it will be performed.
        """
        if self._score is None:
            return self.alignment()
        else:
            return self._score

    def _get_score(self, s1: str, s2: str) -> int:
        if self._subst_mtx is not None:
            try:
                return self._subst_mtx.get(s1, s2)
            except KeyError:
                pass
        if s1 == s2:
            return self.match_score
        elif s1 != s2:
            return self.mismatch_score

    def _build_smith_waterman_matrix(self):
        l1 = len(self.seq1) + 1
        l2 = len(self.seq2) + 1
        self._sw_matrix = np.zeros((l1, l2), dtype=int)
        for i in range(1, l1):
            for j in range(1, l2):
                this_score = self._get_score(self.seq1[i - 1], self.seq2[j - 1])
                if self.seq1[i - 1] == self.seq2[j - 1]:
                    score = self._sw_matrix[i - 1][j - 1] + this_score
                else:
                    score = max(
                        self._sw_matrix[i - 1][j] + self.indel_score,
                        self._sw_matrix[i][j - 1] + self.indel_score,
                        self._sw_matrix[i - 1][j - 1] + this_score,
                    )
                if not self.is_global:
                    score = max(score, 0)
                self._sw_matrix[i][j] = score

    def get_backtrack(
        self,
        seq1_title: str = "seq1",
        seq2_title: str = "seq2",
        alignment_title: str = "aln",
    ) -> Optional[List[str]]:
        """
        Get Smith-Waterman Alignment Backtrack in a human-readable form.
        """
        if self._backtrack is not None:
            return self._backtrack
        if self._sw_matrix is None:
            self._build_smith_waterman_matrix()
        l1 = len(self.seq1)
        l2 = len(self.seq2)

        @functools.lru_cache()
        def location_inside_bound(input_location: Tuple[int, int]) -> bool:
            """Tell whether we're still inside the matrix"""
            return 0 <= input_location[0] <= l1 and 0 <= input_location[1] <= l2

        # warnings.warn(repr(self._sw_matrix))
        best_match_locations = np.where(self._sw_matrix == self.score)
        """All locations with highest scores"""

        retl = set()
        """Return list of strings, one for each best path"""

        for best_match_location in zip(*best_match_locations):
            # Assertion to ensure data correctness
            if len(best_match_location) != 2:
                continue
            # Find a route to origin
            location = best_match_location
            step = []
            while location_inside_bound(location) and sum(location) != 0:
                step.append(location)
                next_score = -math.inf
                next_location = (-1, -1)
                for next_possible_location in [
                    (location[0] - 1, location[1] - 1),
                    (location[0], location[1] - 1),
                    (location[0] - 1, location[1]),
                ]:
                    if (
                        location_inside_bound(next_possible_location)
                        and self._sw_matrix[next_possible_location] > next_score
                    ):
                        next_score = self._sw_matrix[next_possible_location]
                        next_location = next_possible_location
                location = next_location
            # Find a route to maximum point
            location = best_match_location
            step.reverse()
            while location_inside_bound(location) and sum(location) != l1 + l2:
                step.append(location)
                next_score = -math.inf
                next_location = (-1, -1)
                for next_possible_location in [
                    (location[0] + 1, location[1] + 1),
                    (location[0], location[1] + 1),
                    (location[0] + 1, location[1]),
                ]:
                    if (
                        location_inside_bound(next_possible_location)
                        and self._sw_matrix[next_possible_location] > next_score
                    ):
                        next_score = self._sw_matrix[next_possible_location]
                        next_location = next_possible_location
                location = next_location
            step.append((l1, l2))
            # warnings.warn(repr(step))
            # Skip if too few steps
            if len(step) <= 2:
                continue

            # Start compiling strings
            prev_location = (0, 0)
            out_array = []
            for this_location in step:
                if this_location == prev_location:
                    continue
                elif this_location == (prev_location[0] + 1, prev_location[1] + 1):
                    if self.seq1[this_location[0] - 1] == self.seq2[this_location[1] - 1]:
                        out_array.append((self.seq1[this_location[0] - 1], "=", self.seq2[this_location[1] - 1]))
                    else:
                        out_array.append((self.seq1[this_location[0] - 1], "M", self.seq2[this_location[1] - 1]))
                elif this_location == (prev_location[0] + 1, prev_location[1]):
                    out_array.append((self.seq1[this_location[0] - 1], "D", "-"))
                elif this_location == (prev_location[0], prev_location[1] + 1):
                    out_array.append(("-", "I", self.seq2[this_location[1] - 1]))
                else:
                    out_array.append(("-", "?", "-"))
                prev_location = this_location
            rets = (
                ":".join((f">{alignment_title}", seq1_title, "qual", seq2_title, str(np.max(self._sw_matrix[1:, 1:]))))
                + "\n"
                + "\n".join((_seq1 + "" for _seq1 in ["".join(_seq) for _seq in zip(*out_array)]))
            )
            retl.add(rets)
        self._backtrack = list(retl)
        return self._backtrack


def hamming_distance(str1: str, str2: str) -> int:
    """
    Generate hamming distance.

    :raise ValueError: raise this error if the length of input string is not equal.

    >>> hamming_distance("AAAA", "AATA")
    1

    .. versionadded:: 1.0.2
    """
    if len(str1) != len(str2):
        raise ValueError(f"Length of input string 1={str1}, 2={str2} is not equal.")
    return sum(el1 != el2 for el1, el2 in zip(str1, str2))


def editing_distance(str1: str, str2: str) -> int:
    """
    Generate editing distance

    >>> editing_distance("AAAA", "AAAA")
    0
    >>> editing_distance("AAAA", "AATA")
    1
    >>> editing_distance("AAAA", "AAA")
    1
    >>> editing_distance("AAAA", "TAAAA")
    1
    >>> editing_distance("A", "A")
    0
    >>> editing_distance("", "")
    0
    >>> editing_distance("", "A")
    1

    .. versionadded:: 1.0.2
    """
    l1 = len(str1)
    l2 = len(str2)
    if l1 == 0:
        return l2
    elif l2 == 0:
        return l1
    l1 += 1
    l2 += 1
    score_matrix = np.zeros((l1, l2), dtype=int)
    score_matrix[0, :] = range(l2)
    score_matrix[:, 0] = range(l1)
    for i in range(1, l1):
        for j in range(1, l2):
            if str1[i - 1] == str2[j - 1]:
                score_matrix[i][j] = score_matrix[i - 1][j - 1]
            else:
                score_matrix[i][j] = (
                    np.min(
                        (
                            score_matrix[i - 1][j - 1],
                            score_matrix[i][j - 1],
                            score_matrix[i - 1][j],
                        )
                    )
                    + 1
                )
    return score_matrix[l1 - 1][l2 - 1].item()
