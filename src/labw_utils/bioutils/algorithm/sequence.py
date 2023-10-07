"""
``labw_utils.bioutils.algorithm.sequence`` -- Naive sequence algorithms. e.g., complement, reverse or get GC content.

.. versionadded:: 1.0.2
"""
import itertools

from labw_utils.typing_importer import Iterable, List, Optional, Tuple

_comp_trans = str.maketrans("ATCGatcgNnXx", "TAGCtagcNnXx")

AA_SYMBOLS = {
    "A": "Ala",
    "B": "Asx",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
    "Z": "Glx",
    "X": "Any",
    "*": "Stp",
}
"""
TODO: docs

.. versionadded:: 1.0.2
"""

AA_NAMES = {
    "A": "Alanine",
    "B": "Asparagine or Aspartic acid",
    "C": "Cysteine",
    "D": "Aspartic acid",
    "E": "Glutamic acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "Z": "Glutamine or Glutamic acid",
    "X": "Any",
    "*": "Stop",
}
"""
TODO: docs

.. versionadded:: 1.0.2
"""


TRANSL_TABLES_NT = list(map("".join, itertools.product(*["TCAG"] * 3)))
"""
TODO: docs

.. versionadded:: 1.0.2
"""


TRANSL_TABLES = {
    1: {"AA": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "NAME": "The Standard Code"},
    2: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "NAME": "The Vertebrate Mitochondrial Code",
    },
    3: {
        "AA": "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Yeast Mitochondrial Code",
    },
    4: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
    },
    5: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "NAME": "The Invertebrate Mitochondrial Code",
    },
    6: {
        "AA": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Ciliate, Dasycladacean and Hexamita Nuclear Code",
    },
    9: {
        "AA": " FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "NAME": "The Echinoderm and Flatworm Mitochondrial Code",
    },
    10: {"AA": "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "NAME": "The Euplotid Nuclear Code"},
    11: {
        "AA": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Bacterial, Archaeal and Plant Plastid Code",
    },
    12: {
        "AA": "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Alternative Yeast Nuclear Code",
    },
    13: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        "NAME": "The Ascidian Mitochondrial Code",
    },
    14: {
        "AA": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "NAME": "The Alternative Flatworm Mitochondrial Code",
    },
    16: {
        "AA": "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Chlorophycean Mitochondrial Code",
    },
    21: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "NAME": "Trematode Mitochondrial Code",
    },
    22: {
        "AA": "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Scenedesmus obliquus Mitochondrial Code",
    },
    23: {
        "AA": "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Thraustochytrium Mitochondrial Code",
    },
    24: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        "NAME": "Rhabdopleuridae Mitochondrial Code",
    },
    25: {
        "AA": "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Candidate Division SR1 and Gracilibacteria Code",
    },
    26: {
        "AA": "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Pachysolen tannophilus Nuclear Code",
    },
    27: {"AA": "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "NAME": "Karyorelict Nuclear Code"},
    28: {"AA": "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "NAME": "Condylostoma Nuclear Code"},
    29: {"AA": "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "NAME": "Mesodinium Nuclear Code"},
    30: {"AA": "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "NAME": "Peritrich Nuclear Code"},
    31: {
        "AA": "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Blastocrithidia Nuclear Code",
    },
    33: {
        "AA": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        "NAME": "Cephalodiscidae Mitochondrial UAA-Tyr Code",
    },
}
"""
NCBI Translation Table.

This table does NOT care circumstances that a codon may encode either normal AA or STOP.

.. versionadded:: 1.0.2
"""


class MalformedMRNAError(ValueError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    pass


def is_valid_chrname(chr_name: str) -> bool:
    """
    Whether this chrname is a chromosome instead of some alt, ref, patch, etc.

    .. versionadded:: 1.0.2
    """
    if chr_name.startswith("N"):
        if not chr_name.startswith("NC_"):
            return False
    elif chr_name.startswith("chr"):
        if (
            chr_name.startswith("chrUn")
            or chr_name.endswith("_random")
            or chr_name.endswith("_alt")
            or chr_name.endswith("_decoy")
            or chr_name.endswith("_fix")
        ):
            return False
    else:
        if (
            chr_name.startswith("KI")
            or chr_name.startswith("KB")
            or chr_name.startswith("KN")
            or chr_name.startswith("KV")
            or chr_name.startswith("KZ")
            or chr_name.startswith("ML")
            or chr_name.startswith("GL")
            or chr_name.startswith("KQ")
            or chr_name.startswith("CHR_")
            or chr_name.startswith("JH")
        ):
            return False
    return True


def find_orf(seq: str, transl_table: int = 1, init_codon: Optional[List[str]] = None) -> Iterable[Tuple[int, int]]:
    """
    :param seq: The sequence.
    :param transl_table: NCBI Translation Table.
        See <https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes> for more details.
    :param init_codon: Possible initiation codon.
        This is different in different organisms. You should referr to documentations at ``transl_table``.

    .. versionadded:: 1.0.2
    """
    raise NotImplementedError


def translate_cdna(seq: str, transl_table: int = 1) -> str:
    """
    :param seq: The sequence.
    :param transl_table: NCBI Translation Table.
        See <https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes> for more details.

    Translate RNA to AA.
    >>> translate_cdna("TACCGGGTTAATAGGAAACTGACATTTGGAGCCAACACTAGAGGAATCATGAAACTC")
    'YRVNRKLTFGANTRGIMKL'

     .. versionadded:: 1.0.2
    """
    if seq == "":
        return ""
    if len(seq) < 3:
        raise MalformedMRNAError(f"seq ('{seq}') too short")
    if len(seq) % 3 != 0:
        raise MalformedMRNAError(f"Length of seq ('{len(seq)}') should be a multiple of 3")
    seq = seq.upper()
    aa_table = TRANSL_TABLES[transl_table]["AA"]
    return "".join(
        aa_table[TRANSL_TABLES_NT.index(seq[i : i + 3])] if "N" not in seq[i : i + 3] else "X"
        for i in range(0, len(seq), 3)
    )


def complement(seq: str) -> str:
    """
    Get complement of a sequence

    >>> complement("CTGACTGA")
    'GACTGACT'

     .. versionadded:: 1.0.2
    """
    return seq.translate(_comp_trans)


def reverse_complement(seq: str) -> str:
    """
    Get reverse-complement of a sequence

    >>> reverse_complement("CTGACTGA")
    'TCAGTCAG'

     .. versionadded:: 1.0.2
    """
    return complement(seq)[::-1]


def get_gc_percent(seq: str) -> float:
    """
    Get GC content.

    >>> get_gc_percent("AAACG")
    0.4

    Implementation details. Code::

        for base in ("C", "G", "c", "g"):
            gc += seq.count(base)

    is slower than::

        for base in seq:
            if base in ("C", "G", "c", "g"):
                gc += 1

    by 12 percent.

     .. versionadded:: 1.0.2
    """
    if seq is None:
        return None
    if len(seq) == 0:
        return 0
    gc = 0
    for base in seq:
        if base in ("C", "G", "c", "g"):
            gc += 1
    return gc / len(seq)


def decode_phred33(seq: str) -> Iterable[int]:
    """
    Decode phred33 scores (Q-scores) in FASTQ files.

     .. versionadded:: 1.0.2
    """
    for i in seq:
        yield ord(i) - 33
