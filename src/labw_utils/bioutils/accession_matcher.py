"""
``labw_utils.bioutils.accession_matcher`` -- Describe accession information.

This module can identify accession information from its name.
You can see source database and other related information, if possible.

Currently supporting:

- NCBI GenBank Nucleotides (NT) and Whole-Shotgun Sequencing (WGS)
- NCBI RefSeq
- Analysis Set (UCSC style) chromosomal identifier.

>>> import pprint
>>> print(infer_accession_type("chr1"))
Analysis Set Chromosome, with {'NUMBER': '1', 'TYPE': 'NORMAL'}
>>> pprint.pprint(infer_accession_type("CM000663.2").as_dict())
{'details': {'DATABASE': 'NCBI',
             'DTYPE': 'NUCLEOTIDES',
             'TYPE': 'CON division',
             'VERSION': '2'},
 'toplevel': 'NCBI GeneBank Sequence'}
>>> print(infer_accession_type("NC_000001.11"))
NCBI RefSeq, with {'TYPE': 'Chromosome in Reference Assembly', 'VERSION': '11'}
>>> print(infer_accession_type("NC_000001"))
NCBI RefSeq, with {'TYPE': 'Chromosome in Reference Assembly', 'VERSION': None}
>>> pprint.pprint(infer_accession_type("chrUn_JTFH01000804v1_decoy").as_dict())
{'details': {'GENBANK_ACCESSION': 'JTFH01000804',
             'GENBANK_MATCH': {'details': {'DATABASE': 'GenBank',
                                           'DTYPE': 'Whole Genome Shotgun '
                                                    '(WGS)',
                                           'INTERNAL_VERSION': '01',
                                           'TYPE': 'WGS',
                                           'VERSION': None},
                               'toplevel': 'NCBI GeneBank Sequence'},
             'TYPE': 'DECOY',
             'VERSION': '1'},
 'toplevel': 'Analysis Set Unplaced Scaffold'}
>>> pprint.pprint(infer_accession_type("chr19_KI270885v1_alt").as_dict())
{'details': {'GENBANK_ACCESSION': 'KI270885',
             'GENBANK_MATCH': {'details': {'DATABASE': 'NCBI',
                                           'DTYPE': 'NUCLEOTIDES',
                                           'TYPE': 'CON division',
                                           'VERSION': None},
                               'toplevel': 'NCBI GeneBank Sequence'},
             'PLACED_CHROMOSOME_NUMBER': '19',
             'TYPE': 'Alternate Loci',
             'VERSION': '1'},
 'toplevel': 'Analysis Set Placed Scaffold'}
>>> pprint.pprint(infer_accession_type("chr22_KI270733v1_random").as_dict())
{'details': {'GENBANK_ACCESSION': 'KI270733',
             'GENBANK_MATCH': {'details': {'DATABASE': 'NCBI',
                                           'DTYPE': 'NUCLEOTIDES',
                                           'TYPE': 'CON division',
                                           'VERSION': None},
                               'toplevel': 'NCBI GeneBank Sequence'},
             'PLACED_CHROMOSOME_NUMBER': '22',
             'TYPE': 'Unlocalized Scaffolds',
             'VERSION': '1'},
 'toplevel': 'Analysis Set Placed Scaffold'}

.. warning::
    Is not finished.

.. versionadded:: 1.0.2
"""

from __future__ import annotations

__all__ = ("infer_accession_type", "AccessionMatchResult")

import re
from abc import abstractmethod, ABC

from labw_utils.devutils.decorators import create_class_init_doc_from_property, supress_inherited_doc
from labw_utils.typing_importer import Union, Optional, Type, Final, Any, Callable, Mapping, Iterable, List


@create_class_init_doc_from_property()
class AccessionMatchResult:
    """
    Accession match result.

    .. versionadded:: 1.0.2
    """

    _toplevel: str
    _details: Mapping[str, AccessionMatchResult]

    def __init__(self, toplevel: str, details: Mapping[str, Union[str, AccessionMatchResult]]):
        self._toplevel = toplevel
        self._details = details

    def as_dict(self) -> Mapping[str, Union[str, Mapping[str, str]]]:
        """Convert to :py:class:`dict`"""
        return {
            "toplevel": self._toplevel,
            "details": {k: v.as_dict() if isinstance(v, AccessionMatchResult) else v for k, v in self._details.items()},
        }

    def __repr__(self) -> str:
        return f"{self._toplevel}, with {self._details}"

    @property
    def toplevel(self) -> str:
        """The toplevel identifier"""
        return self._toplevel

    @property
    def details(self) -> Mapping[str, Union[str, AccessionMatchResult]]:
        """
        Other details. May contain information to match results of sub-accessions.
        """
        return self._details


class AccessionMatcherRuleType(ABC):
    """
    Abstract rules.

    .. versionadded:: 1.0.2
    """

    @abstractmethod
    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        """
        To be overridden by developers.
        """
        raise NotImplementedError


@supress_inherited_doc(modify_overwritten=True)
class ChainAccessionMatcherRuleType(AccessionMatcherRuleType):
    """
    A rule matcher that would try all rules in pre-defined order and report result of the first match.

    .. versionadded:: 1.0.2
    """

    _rule_chain: List[Type[AccessionMatcherRuleType]]

    def __init__(self):
        """
        :raises TypeError: If ``_rule_chain`` class variable was not set.
        """
        if getattr(self, "_rule_chain", None) is None:
            raise TypeError

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        for rule in self._rule_chain:
            match_result = rule().match(accession)
            if match_result is not None:
                return match_result
        return None


@supress_inherited_doc(modify_overwritten=True)
class MasterEnsembleIDMatcher(AccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _regex = re.compile(r"^ENS([A-Z_]+)?(G|T|E|P|R|FM|GT)\d+(\.\d+)?$")

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        if not accession.startswith("ENS"):
            return None
        match_result = MasterEnsembleIDMatcher._regex.match(accession)
        if match_result is None:
            return AccessionMatchResult("Ensemble ID", {})
        groups = match_result.groups()
        details_dict = {}
        if groups[0] is None:
            details_dict["ORGANISM"] = "N/A"
        else:
            details_dict["ORGANISM"] = groups[0]
        details_dict["FEATURE"] = {
            "G": "Gene",
            "T": "Transcript",
            "E": "Exon",
            "FM": "Protein Family",
            "GT": "Gene Tree",
            "P": "Protein",
            "R": "Regulatory Feature",
        }.get(groups[1], "UNKNOWN")
        if groups[2] is None:
            details_dict["VERSION"] = "N/A"
        else:
            details_dict["VERSION"] = groups[2][1:]
        return AccessionMatchResult("Ensemble ID", details_dict)


@supress_inherited_doc(modify_overwritten=True)
class NCBIRefSeqMatcher(AccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _regex_dict = {
        "Protein-Coding Transcript": re.compile(r"^NM_\d+(\.\d+)?$"),
        "Non-Protein-Coding Transcript": re.compile(r"^NR_\d+(\.\d+)?$"),
        "Genomic Contig or scaffold, clone-based or WGS": re.compile(r"^NT_\d+(\.\d+)?$"),
        "Genomic Contig or scaffold, primarily WGS": re.compile(r"^NW_\d+(\.\d+)?$"),
        "Genomic Complete genomes and unfinished WGS data": re.compile(r"^NZ_\d+(\.\d+)?$"),
        "Predicted Protein-Coding Transcript": re.compile(r"^XM_\d+(\.\d+)?$"),
        "Predicted Non-Protein-Coding Transcript": re.compile(r"^XR_\d+(\.\d+)?$"),
        "Chromosome in Reference Assembly": re.compile(r"^NC_\d+(\.\d+)?$"),
        "Chromosome in Alternate Assembly": re.compile(r"^AC_\d+(\.\d+)?$"),
        "Incomplete Genomic Region": re.compile(r"^NG_\d+(\.\d+)?$"),
        "Protein Annotated on AC_ alternate assembly": re.compile(r"^AP_\d+(\.\d+)?$"),
        "Protein Associated with an NM_ or NC_ accession": re.compile(r"^NP_\d+(\.\d+)?$"),
        "Protein Annotated on genomic molecules without an instantiated transcript record": re.compile(
            r"^YP_\d+(\.\d+)?$"
        ),
        "Protein Predicted model, associated with an XM_ accession": re.compile(r"^XP_\d+(\.\d+)?$"),
        "Protein	Non-redundant across multiple strains and species": re.compile(r"^WP_\d+(\.\d+)?$"),
    }

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        for (
            k,
            v,
        ) in NCBIRefSeqMatcher._regex_dict.items():
            match_result = v.match(accession)
            if match_result is not None:
                version = match_result.groups()[0]
                if version is not None:
                    version = version.strip(".")
                return AccessionMatchResult(toplevel="NCBI RefSeq", details={"TYPE": k, "VERSION": version})
        return None


@supress_inherited_doc(modify_overwritten=True)
class AnalysisSetChromosomeHLA(AccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _regex = re.compile(r"^HLA-(.+)$")

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        match_result = self._regex.match(accession)
        if match_result is None:
            return None
        groups = match_result.groups()
        hla_name = groups[0]
        return AccessionMatchResult(toplevel="Analysis Set HLA Sequence", details={"HLA_NAME": hla_name})


@supress_inherited_doc(modify_overwritten=True)
class NCBIGenBankAccessionMatchingEngine:
    """
    TODO: add protein support

    .. versionadded:: 1.0.2
    """

    _regex_dict: Optional[Mapping[re.Pattern, Callable[[str], Mapping[str, Any]]]]

    def __init__(self):
        self._regex_dict = None

    @staticmethod
    def generate_ncbi_genbank_nt_regex(
        prefixes: Iterable[str], kwargs: Mapping[str, str]
    ) -> Mapping[re.Pattern, Callable[[str], Mapping[str, Any]]]:
        def generate_ncbi_genbank_nt_regex_from_single_prefix(
            _prefix: str,
        ) -> Mapping[re.Pattern, Callable[[str], Mapping[str, Any]]]:
            if len(_prefix) == 1:
                regex = r"^" + _prefix + r"\d{5}" + r"(\.\d+)?" + r"$"
            elif len(_prefix) == 2:
                regex = r"^" + _prefix + r"\d{6,8}" + r"(\.\d+)?" + r"$"
            else:
                raise ValueError(_prefix)
            regex = re.compile(regex)

            def retf(accession: str) -> Mapping[str, str]:
                match_result = regex.match(accession)
                if match_result is not None:
                    version = match_result.groups()[0]
                    if version is not None:
                        version = version.strip(".")
                else:
                    raise ValueError
                return {**kwargs, "DTYPE": "NUCLEOTIDES", "VERSION": version}

            return {regex: retf}

        retd = {}
        for prefix in prefixes:
            retd.update(generate_ncbi_genbank_nt_regex_from_single_prefix(prefix))
        return retd

    @staticmethod
    def generate_ncbi_genbank_wgs_regex(
        prefix_regex: str, prefix_len: int, kwargs: Mapping[str, Any]
    ) -> Mapping[re.Pattern, Callable[[str], Mapping[str, Any]]]:
        a_z_regex = r"[A-Z]" + r"{" + str(prefix_len - 1) + r"}"
        if prefix_len == 4:
            regex = r"^" + prefix_regex + a_z_regex + r"(\d{2})\d{6}(\.\d+)?$"
        elif prefix_len == 6:
            regex = r"^" + prefix_regex + a_z_regex + r"(\d{2})\d{7}(\.\d+)?$"
        else:
            raise ValueError(prefix_len)
        regex = re.compile(regex)

        def retf(accession: str) -> Mapping[str, str]:
            match_result = regex.match(accession)
            if match_result is not None:
                internal_version = match_result.groups()[0]
                version = match_result.groups()[1]
                if version is not None:
                    version = version.strip(".")
            else:
                raise ValueError
            return {
                **kwargs,
                "DTYPE": "Whole Genome Shotgun (WGS)",
                "INTERNAL_VERSION": internal_version,
                "VERSION": version,
            }

        return {regex: retf}

    @property
    def regex_dict(self) -> Mapping[re.Pattern, Callable[[str], Mapping[str, Any]]]:
        if self._regex_dict is None:
            self._regex_dict = {
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("BA", "DF", "DG", "LD"), {"DATABASE": "DDBJ", "TYPE": "CON division"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AN",), {"DATABASE": "EMBL", "TYPE": "CON division"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "CH",
                        "CM",
                        "DS",
                        "EM",
                        "EN",
                        "EP",
                        "EQ",
                        "FA",
                        "GG",
                        "GL",
                        "JH",
                        "KB",
                        "KD",
                        "KE",
                        "KI",
                        "KK",
                        "KL",
                        "KN",
                        "KQ",
                        "KV",
                        "KZ",
                        "ML",
                        "MU",
                    ),
                    {"DATABASE": "NCBI", "TYPE": "CON division"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "C",
                        "AT",
                        "AU",
                        "AV",
                        "BB",
                        "BJ",
                        "BP",
                        "BW",
                        "BY",
                        "CI",
                        "CJ",
                        "DA",
                        "DB",
                        "DC",
                        "DK",
                        "FS",
                        "FY",
                        "HX",
                        "HY",
                        "LU",
                        "OH",
                    ),
                    {"DATABASE": "DDBJ", "TYPE": "EST"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("F",), {"DATABASE": "EMBL", "TYPE": "EST"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "H",
                        "N",
                        "T",
                        "R",
                        "W",
                        "AA",
                        "AI",
                        "AW",
                        "BE",
                        "BF",
                        "BG",
                        "BI",
                        "BM",
                        "BQ",
                        "BU",
                        "CA",
                        "CB",
                        "CD",
                        "CF",
                        "CK",
                        "CN",
                        "CO",
                        "CV",
                        "CX",
                        "DN",
                        "DR",
                        "DT",
                        "DV",
                        "DW",
                        "DY",
                        "EB",
                        "EC",
                        "EE",
                        "EG",
                        "EH",
                        "EL",
                        "ES",
                        "EV",
                        "EW",
                        "EX",
                        "EY",
                        "FC",
                        "FD",
                        "FE",
                        "FF",
                        "FG",
                        "FK",
                        "FL",
                        "GD",
                        "GE",
                        "GH",
                        "GO",
                        "GR",
                        "GT",
                        "GW",
                        "HO",
                        "HS",
                        "JG",
                        "JK",
                        "JZ",
                    ),
                    {"DATABASE": "GenBank", "TYPE": "EST"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("D", "AB", "LC"), {"DATABASE": "DDBJ", "TYPE": "Direct submissions"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "V",
                        "X",
                        "Y",
                        "Z",
                        "AJ",
                        "AM",
                        "FM",
                        "FN",
                        "HE",
                        "HF",
                        "HG",
                        "FO",
                        "LK",
                        "LL",
                        "LM",
                        "LN",
                        "LO",
                        "LR",
                        "LS",
                        "LT",
                        "OA",
                        "OB",
                        "OC",
                        "OD",
                        "OE",
                        "OU",
                        "OV",
                        "OW",
                        "OX",
                        "OY",
                        "OZ",
                    ),
                    {"DATABASE": "EMBL", "TYPE": "Direct submissions"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "U",
                        "AF",
                        "AY",
                        "DQ",
                        "EF",
                        "EU",
                        "FJ",
                        "GQ",
                        "GU",
                        "HM",
                        "HQ",
                        "JF",
                        "JN",
                        "JQ",
                        "JX",
                        "KC",
                        "KF",
                        "KJ",
                        "KM",
                        "KP",
                        "KR",
                        "KT",
                        "KU",
                        "KX",
                        "KY",
                        "MF",
                        "MG",
                        "MH",
                        "MK",
                        "MN",
                        "MT",
                        "MW",
                        "MZ",
                        "OK",
                        "OL",
                        "OM",
                        "ON",
                        "OP",
                        "OQ",
                    ),
                    {"DATABASE": "GenBank", "TYPE": "Direct submissions"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AP", "BS"), {"DATABASE": "DDBJ", "TYPE": "Genome project data"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AL", "BX", "CR", "CT", "CU", "FP", "FQ", "FR"),
                    {"DATABASE": "EMBL", "TYPE": "Genome project data"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AE", "CP", "CY"), {"DATABASE": "GenBank", "TYPE": "Genome project data"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AG", "DE", "DH", "FT", "GA", "LB"), {"DATABASE": "DDBJ", "TYPE": "GSS"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "B",
                        "AQ",
                        "AZ",
                        "BH",
                        "BZ",
                        "CC",
                        "CE",
                        "CG",
                        "CL",
                        "CW",
                        "CZ",
                        "DU",
                        "DX",
                        "ED",
                        "EI",
                        "EJ",
                        "EK",
                        "ER",
                        "ET",
                        "FH",
                        "FI",
                        "GS",
                        "HN",
                        "HR",
                        "JJ",
                        "JM",
                        "JS",
                        "JY",
                        "KG",
                        "KO",
                        "KS",
                        "MJ",
                    ),
                    {"DATABASE": "GenBank", "TYPE": "GSS"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AK"), {"DATABASE": "DDBJ", "TYPE": "cDNA projects"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AC", "DP"), {"DATABASE": "GenBank", "TYPE": "HTGS"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "E",
                        "BD",
                        "DD",
                        "DI",
                        "DJ",
                        "DL",
                        "DM",
                        "FU",
                        "FV",
                        "FW",
                        "FZ",
                        "GB",
                        "HV",
                        "HW",
                        "HZ",
                        "LF",
                        "LG",
                        "LV",
                        "LX",
                        "LY",
                        "LZ",
                        "MA",
                        "MB",
                        "MC",
                        "MD",
                        "ME",
                        "OF",
                        "OG",
                        "OI",
                        "OJ",
                        "PA",
                        "PB",
                        "PC",
                        "PD",
                        "PE",
                    ),
                    {"DATABASE": "DDBJ", "TYPE": "Patents"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "A",
                        "AX",
                        "CQ",
                        "CS",
                        "FB",
                        "GM",
                        "GN",
                        "HA",
                        "HB",
                        "HC",
                        "HD",
                        "HH",
                        "HI",
                        "JA",
                        "JB",
                        "JC",
                        "JD",
                        "JE",
                        "LP",
                        "LQ",
                        "MP",
                        "MQ",
                        "MR",
                        "MS",
                    ),
                    {"DATABASE": "EMBL", "TYPE": "Patents (nucleotide only)"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "I",
                        "AR",
                        "DZ",
                        "EA",
                        "GC",
                        "GP",
                        "GV",
                        "GX",
                        "GY",
                        "GZ",
                        "HJ",
                        "HK",
                        "HL",
                        "KH",
                        "MI",
                        "MM",
                        "MO",
                        "MV",
                        "MX",
                        "MY",
                        "OO",
                    ),
                    {"DATABASE": "GenBank", "TYPE": "Patents (nucleotide)"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "G",
                        "BV",
                        "GF",
                    ),
                    {"DATABASE": "GenBank", "TYPE": "STS"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("BR",), {"DATABASE": "DDBJ", "TYPE": "TPA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("BN",), {"DATABASE": "EMBL", "TYPE": "TPA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("BK",), {"DATABASE": "GenBank", "TYPE": "TPA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "HT",
                        "HU",
                    ),
                    {"DATABASE": "DDBJ", "TYPE": "TPA CON division"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "BL",
                        "GJ",
                        "GK",
                    ),
                    {"DATABASE": "GenBank", "TYPE": "TPA CON division"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "EZ",
                        "HP",
                        "JI",
                        "JL",
                        "JO",
                        "JP",
                        "JR",
                        "JT",
                        "JU",
                        "JV",
                        "JW",
                        "KA",
                    ),
                    {"DATABASE": "GenBank", "TYPE": "TSA"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "FX",
                        "LA",
                        "LE",
                        "LH",
                        "LI",
                        "LJ",
                    ),
                    {"DATABASE": "DDBJ", "TYPE": "TSA"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("S",), {"DATABASE": "GenBank", "TYPE": "From journal scanning"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AD",), {"DATABASE": "GenBank", "TYPE": "From GSDB"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AH",), {"DATABASE": "GenBank", "TYPE": "Segmented set header"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("AS",), {"DATABASE": "GenBank", "TYPE": "Other - not currently being used"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("BC",), {"DATABASE": "GenBank", "TYPE": "MGC project"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    ("BT",), {"DATABASE": "GenBank", "TYPE": "FLI-cDNA projects"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_nt_regex(
                    (
                        "J",
                        "K",
                        "L",
                        "M",
                    ),
                    {"DATABASE": "GenBank", "TYPE": "From GSDB direct submissions"},
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"[AJLMMPQRSVWX]", 4, {"DATABASE": "GenBank", "TYPE": "WGS"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"[AJ]", 6, {"DATABASE": "GenBank", "TYPE": "WGS"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"B", 4, {"DATABASE": "DDBJ", "TYPE": "WGS"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"B", 6, {"DATABASE": "DDBJ", "TYPE": "WGS"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"[CFOU]", 4, {"DATABASE": "EMBL", "TYPE": "WGS"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"C", 6, {"DATABASE": "EMBL", "TYPE": "WGS"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"D", 4, {"DATABASE": "GenBank", "TYPE": "WGS or TSA TPA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"D", 6, {"DATABASE": "GenBank", "TYPE": "WGS TPA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"E", 4, {"DATABASE": "DDBJ", "TYPE": "WGS TPA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"G", 4, {"DATABASE": "GenBank", "TYPE": "TSA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"H", 4, {"DATABASE": "EMBL", "TYPE": "TSA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"I", 4, {"DATABASE": "DDBJ", "TYPE": "TSA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"T", 4, {"DATABASE": "DDBJ", "TYPE": "Targeted Gene Projects"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"K", 4, {"DATABASE": "GenBank", "TYPE": "Targeted Gene Projects (TLS)"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"Y", 4, {"DATABASE": "DDBJ", "TYPE": "TSA TPA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"Z", 4, {"DATABASE": "DDBJ", "TYPE": "Targeted Gene Projects TPA"}
                ),
                **NCBIGenBankAccessionMatchingEngine.generate_ncbi_genbank_wgs_regex(
                    r"A", 6, {"DATABASE": "DDBJ", "TYPE": "MGA"}
                ),
            }
        return self._regex_dict


ncbi_genbank_accession_matching_engine = NCBIGenBankAccessionMatchingEngine()


@supress_inherited_doc(modify_overwritten=True)
class NCBIGenBankAccessionMatcher(AccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    def __init__(self):
        self._regex_dict = ncbi_genbank_accession_matching_engine.regex_dict

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        for k, v in self._regex_dict.items():
            if k.match(accession) is not None:
                return AccessionMatchResult(toplevel="NCBI GeneBank Sequence", details=v(accession))


@supress_inherited_doc(modify_overwritten=True)
class AnalysisSetChromosomePlacedGenbankMatcher(AccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _regex = re.compile(r"^chr(\d+|X|Y|M)_(.+)v(\d+)?(_alt|_random|_fix)?$")

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        match_result = self._regex.match(accession)
        if match_result is None:
            return None
        groups = match_result.groups()
        genbank_accession = groups[1]
        return AccessionMatchResult(
            toplevel="Analysis Set Placed Scaffold",
            details={
                "PLACED_CHROMOSOME_NUMBER": groups[0],
                "GENBANK_ACCESSION": genbank_accession,
                "GENBANK_MATCH": NCBIGenBankAccessionMatcher().match(genbank_accession),
                "VERSION": groups[2],
                "TYPE": {"_alt": "Alternate Loci", "_random": "Unlocalized Scaffolds", "_fix": "Fix Patches"}.get(
                    groups[3], "N/A"
                ),
            },
        )


@supress_inherited_doc(modify_overwritten=True)
class AnalysisSetChromosomeUnplacedGenbankMatcher(AccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _regex = re.compile(r"^chrUn_(.+)v(\d+)(_decoy)?$")

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        match_result = self._regex.match(accession)
        if match_result is None:
            return None
        groups = match_result.groups()
        genbank_accession = groups[0]
        return AccessionMatchResult(
            toplevel="Analysis Set Unplaced Scaffold",
            details={
                "GENBANK_ACCESSION": genbank_accession,
                "GENBANK_MATCH": NCBIGenBankAccessionMatcher().match(genbank_accession),
                "VERSION": groups[1],
                "TYPE": "DECOY" if groups[2] is not None else "MISC",
            },
        )


@supress_inherited_doc(modify_overwritten=True)
class AnalysisSetChromosomeMatcher(AccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _regex = re.compile(r"^chr(\d+|X|Y|M|EBV)$")

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        match_result = self._regex.match(accession)
        if match_result is None:
            return None
        groups = match_result.groups()
        return AccessionMatchResult(
            toplevel="Analysis Set Chromosome",
            details={
                "NUMBER": groups[0],
                "TYPE": {"X": "Sex", "Y": "Sex", "M": "Mitochondria", "EBV": "Epstein-Barr Virus"}.get(
                    groups[0], "NORMAL"
                ),
            },
        )


@supress_inherited_doc(modify_overwritten=True)
class AnalysisSetContigMatcher(ChainAccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _rule_chain: Final[list[AccessionMatcherRuleType]] = [
        AnalysisSetChromosomeMatcher,
        AnalysisSetChromosomeUnplacedGenbankMatcher,
        AnalysisSetChromosomePlacedGenbankMatcher,
        AnalysisSetChromosomeHLA,
    ]


@supress_inherited_doc(modify_overwritten=True)
class IntegerMatcher(AccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _int_regex = re.compile(r"^\d+$")
    _roman_int_regex = re.compile(r"^(?=[MDCLXVI])M*(C[MD]|D?C{0,3})(X[CL]|L?X{0,3})(I[XV]|V?I{0,3})$")

    # See: <https://www.oreilly.com/library/view/regular-expressions-cookbook/9780596802837/ch06s09.html>

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        if self._int_regex.match(accession) is not None:
            return AccessionMatchResult(toplevel="Integer", details={})
        if self._roman_int_regex.match(accession) is not None:
            return AccessionMatchResult(toplevel="Integer (Roman)", details={})
        return None


@supress_inherited_doc(modify_overwritten=True)
class MasterAccessionMatcher(ChainAccessionMatcherRuleType):
    """

    .. versionadded:: 1.0.2
    """

    _rule_chain: Final[list[AccessionMatcherRuleType]] = [
        MasterEnsembleIDMatcher,
        NCBIRefSeqMatcher,
        AnalysisSetContigMatcher,
        NCBIGenBankAccessionMatcher,
        IntegerMatcher,
    ]


def infer_accession_type(accession: str) -> Optional[AccessionMatchResult]:
    """
    :param accession: Input accession.
    :return: Identified information. :py:obj:`None` if failed.

    .. versionadded:: 1.0.2
    """
    mam = MasterAccessionMatcher()
    return mam.match(accession)
