"""
Naive sequence algorithms. e.g., complement, reverse or get GC content.
"""
import itertools
import re
from typing import Iterable, List, Optional, Tuple, Dict

_comp_trans = str.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')

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
    "*": "Stp"
}

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
    "*": "Stop"
}

TRANSL_TABLES_NT = list(map("".join, itertools.product(*["TCAG"] * 3)))

TRANSL_TABLES = {
    1: {
        "AA": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Standard Code"
    },
    2: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "NAME": "The Vertebrate Mitochondrial Code"
    },
    3: {
        "AA": "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Yeast Mitochondrial Code"
    },
    4: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code"
    },
    5: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "NAME": "The Invertebrate Mitochondrial Code"
    },
    6: {
        "AA": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Ciliate, Dasycladacean and Hexamita Nuclear Code"
    },
    9: {
        "AA": " FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "NAME": "The Echinoderm and Flatworm Mitochondrial Code"
    },
    10: {
        "AA": "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Euplotid Nuclear Code"
    },
    11: {
        "AA": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Bacterial, Archaeal and Plant Plastid Code"
    },
    12: {
        "AA": "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "The Alternative Yeast Nuclear Code"
    },
    13: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        "NAME": "The Ascidian Mitochondrial Code"
    },
    14: {
        "AA": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "NAME": "The Alternative Flatworm Mitochondrial Code"
    },
    16: {
        "AA": "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Chlorophycean Mitochondrial Code"
    },
    21: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "NAME": "Trematode Mitochondrial Code"
    },
    22: {
        "AA": "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Scenedesmus obliquus Mitochondrial Code"
    },
    23: {
        "AA": "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Thraustochytrium Mitochondrial Code"
    },
    24: {
        "AA": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        "NAME": "Rhabdopleuridae Mitochondrial Code"
    },
    25: {
        "AA": "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Candidate Division SR1 and Gracilibacteria Code"
    },
    26: {
        "AA": "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Pachysolen tannophilus Nuclear Code"
    },
    27: {
        "AA": "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Karyorelict Nuclear Code"
    },
    28: {
        "AA": "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Condylostoma Nuclear Code"
    },
    29: {
        "AA": "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Mesodinium Nuclear Code"
    },
    30: {
        "AA": "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Peritrich Nuclear Code"
    },
    31: {
        "AA": "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "NAME": "Blastocrithidia Nuclear Code"
    },
    33: {
        "AA": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        "NAME": "Cephalodiscidae Mitochondrial UAA-Tyr Code"
    }
}
"""
NCBI Translation Table.

This table does NOT care circumstances that a codon may encode either normal AA or STOP.
"""

ACCESSION_MATCHER: Dict[str, re.Pattern] = {
    "Ensemble Gene ID": re.compile(r"^ENS([A-Z_]+)?G\d{8}$"),
    "Ensemble Gene ID with Version": re.compile(r"^ENS([A-Z_]+)?G\d{8}\.\d+$"),
    "Ensemble Transcript ID": re.compile(r"^ENS([A-Z_]+)?T\d{8}$"),
    "Ensemble Transcript ID with Version": re.compile(r"^ENS([A-Z_]+)?T\d{8}\.\d+$"),
    "Ensemble Exon ID": re.compile(r"^ENS([A-Z_]+)?E\d{8}$"),
    "Ensemble Exon ID with Version": re.compile(r"^ENS([A-Z_]+)?E\d{8}\.\d+$"),
    "Ensemble Protein Family ID": re.compile(r"^ENS([A-Z_]+)?FM\d{8}$"),
    "Ensemble Protein Family ID with Version": re.compile(r"^ENS([A-Z_]+)?FM\d{8}\.\d+$"),
    "Ensemble Gene Tree ID": re.compile(r"^ENS([A-Z_]+)?GT\d{8}$"),
    "Ensemble Gene Tree ID with Version": re.compile(r"^ENS([A-Z_]+)?GT\d{8}\.\d+$"),
    "Ensemble Protein ID": re.compile(r"^ENS([A-Z_]+)?P\d{8}$"),
    "Ensemble Protein ID with Version": re.compile(r"^ENS([A-Z_]+)?P\d{8}\.\d+$"),
    "Ensemble Regulatory Feature ID": re.compile(r"^ENS([A-Z_]+)?R\d{8}$"),
    "Ensemble Regulatory Feature ID with Version": re.compile(r"^ENS([A-Z_]+)?R\d{8}\.\d+$"),
    "Ensemble Unknown": re.compile(r"^ENS.+$"),
    "RefSeq Protein-Coding Transcript": re.compile(r"^NM_.+$"),
    "RefSeq Non-Protein-Coding Transcript": re.compile(r"^NR_.+$"),
    "RefSeq Predicted Protein-Coding Transcript": re.compile(r"^XM_c.+$"),
    "RefSeq Predicted Non-Protein-Coding Transcript": re.compile(r"^XR_c.+$"),
    "RefSeq Chromosome in Reference Assembly": re.compile(r"^NC_.+$"),
    "RefSeq Chromosome in Alternate Assembly": re.compile(r"^AC_.+$"),
    "RefSeq Incomplete Genomic Region": re.compile(r"^NG_.+$"),
    "Analysis Set Chromosome": re.compile(r"^chr(([0-9]+)|X|Y|M)$"),
    "Analysis Set Chromosome, EBV": re.compile(r"^chrEBV$"),
    "Analysis Set Chromosome, Unplaced Decoy Sequences": re.compile(r"^chrUn_.+_decoy$"),
    "Analysis Set Chromosome, HLA Sequences": re.compile(r"^HLA-.*$"),
    "Analysis Set Chromosome, Unplaced Scaffolds": re.compile(r"^chrUn_.+$"),
    "Analysis Set Chromosome, Alternate Loci": re.compile(r"^chr(([0-9]+)|X|Y|M)_.+_alt$"),
    "Analysis Set Chromosome, Unlocalized Scaffolds": re.compile(r"^chr(([0-9]+)|X|Y|M)_.+_random$"),
    "Analysis Set Chromosome, Fix Patches": re.compile(r"^chr(([0-9]+)|X|Y|M)_.+_fix$"),
    "Analysis Set Chromosome, Unknown": re.compile(r"^chr.+"),
    "Unknown": re.compile(r".*")
}


# NT_	Genomic	Contig or scaffold, clone-based or WGSa
# NW_	Genomic	Contig or scaffold, primarily WGSa
# NZ_b	Genomic	Complete genomes and unfinished WGS data
# AP_	Protein	Annotated on AC_ alternate assembly
# NP_	Protein	Associated with an NM_ or NC_ accession
# YP_c	Protein	Annotated on genomic molecules without an instantiated
# transcript record
# XP_c	Protein	Predicted model, associated with an XM_ accession
# WP_	Protein	Non-redundant across multiple strains and species

class MalformedMRNAError(ValueError):
    pass


def infer_accession_type(accession: str) -> str:
    for k, v in ACCESSION_MATCHER.items():
        if v.match(accession):
            return k


def is_valid_chrname(chr_name: str) -> bool:
    """Whether this chrname is a chromosome instead of some alt, ref, patch, etc."""
    if chr_name.startswith("N"):
        if not chr_name.startswith("NC_"):
            return False
    elif chr_name.startswith("chr"):
        if (
                chr_name.startswith("chrUn") or
                chr_name.endswith("_random") or
                chr_name.endswith("_alt") or
                chr_name.endswith("_decoy") or
                chr_name.endswith("_fix")
        ):
            return False
    else:
        if (
                chr_name.startswith("KI") or
                chr_name.startswith("KB") or
                chr_name.startswith("KN") or
                chr_name.startswith("KV") or
                chr_name.startswith("KZ") or
                chr_name.startswith("ML") or
                chr_name.startswith("GL") or
                chr_name.startswith("KQ") or
                chr_name.startswith("CHR_") or
                chr_name.startswith("JH")
        ):
            return False
    return True


def find_orf(
        seq: str,
        transl_table: int = 1,
        init_codon: Optional[List[str]] = None
) -> Iterable[Tuple[int, int]]:
    """
    :param seq: The sequence.
    :param transl_table: NCBI Translation Table.
        See <https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes> for more details.
    :param init_codon: Possible initiation codon.
        This is different in different organisms. You should referr to documentations at ``transl_table``.

    """
    raise NotImplementedError


def translate_cdna(
        seq: str,
        transl_table: int = 1
) -> str:
    """
    :param seq: The sequence.
    :param transl_table: NCBI Translation Table.
        See <https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes> for more details.

    Translate RNA to AA.
    >>> translate_cdna("TACCGGGTTAATAGGAAACTGACATTTGGAGCCAACACTAGAGGAATCATGAAACTC")
    'YRVNRKLTFGANTRGIMKL'
    """
    if len(seq) < 3:
        raise MalformedMRNAError(f"seq ('{seq}') too short")
    if len(seq) % 3 != 0:
        raise MalformedMRNAError(f"Length of seq ('{len(seq)}') should be a multiple of 3")
    seq = seq.upper()
    aa_table = TRANSL_TABLES[transl_table]["AA"]
    return "".join(
        aa_table[TRANSL_TABLES_NT.index(seq[i:i + 3])]
        for i in range(0, len(seq), 3)
    )


def complement(seq: str) -> str:
    """
    Get complement of a sequence

    >>> complement("CTGACTGA")
    'GACTGACT'
    """
    return seq.translate(_comp_trans)


def reverse_complement(seq: str) -> str:
    """
    Get reverse-complement of a sequence

    >>> reverse_complement("CTGACTGA")
    'TCAGTCAG'
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
    """
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
    """
    for i in seq:
        yield ord(i) - 33
