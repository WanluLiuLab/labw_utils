from __future__ import annotations

import re
from typing import List, Final, Optional, Mapping, Iterable

from labw_utils.bioutils.accession_matcher import ChainAccessionMatcherRuleType, AccessionMatcherRuleType, \
    AccessionMatchResult


def generate_ncbi_genbank_nt_regex(
        prefixes: Iterable[str],
        kwargs: Mapping[str, str]
) -> Mapping[re.Pattern, Mapping[str, str]]:
    retd = {}
    for prefix in prefixes:
        if len(prefix) == 1:
            regex = r"^" + prefix + r"\d{5}" + r"$"
        elif len(prefix) == 2:
            regex = r"^" + prefix + r"\d{6,8}" + r"$"
        else:
            raise ValueError(prefix)
        retd[re.compile(regex)] = kwargs
    return retd


class MasterEnsembleIDMatcher(AccessionMatcherRuleType):
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
            "R": "Regulatory Feature"
        }.get(groups[1], "UNKNOWN")
        if groups[2] is None:
            details_dict["VERSION"] = "N/A"
        else:
            details_dict["VERSION"] = groups[2][1:]
        return AccessionMatchResult("Ensemble ID", details_dict)


class NCBIRefSeqMatcher(AccessionMatcherRuleType):
    _regex_dict = {
        "Protein-Coding Transcript": re.compile(r"^NM_.+$"),
        "Non-Protein-Coding Transcript": re.compile(r"^NR_.+$"),
        "Predicted Protein-Coding Transcript": re.compile(r"^XM_c.+$"),
        "Predicted Non-Protein-Coding Transcript": re.compile(r"^XR_c.+$"),
        "Chromosome in Reference Assembly": re.compile(r"^NC_.+$"),
        "Chromosome in Alternate Assembly": re.compile(r"^AC_.+$"),
        "Incomplete Genomic Region": re.compile(r"^NG_.+$"),
    }

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        for k, v, in NCBIRefSeqMatcher._regex_dict.items():
            if v.match(accession) is not None:
                return AccessionMatchResult(
                    toplevel="NCBI RefSeq",
                    details={"TYPE": k}
                )
        return None


class AnalysisSetChromosomeHLA(AccessionMatcherRuleType):
    _regex = re.compile(r"^HLA-(.+)$")

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        match_result = self._regex.match(accession)
        if match_result is None:
            return None
        groups = match_result.groups()
        hla_name = groups[0]
        return AccessionMatchResult(
            toplevel="Analysis Set HLA Sequence",
            details={
                "HLA_NAME": hla_name
            }
        )


class NCBIGenBankAccessionMatcher(AccessionMatcherRuleType):
    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        for k, v in self._regex_dict.items():
            if k.match(accession):
                return AccessionMatchResult(
                    toplevel="NCBI GeneBank Sequence",
                    details=v
                )

    _regex_dict: Mapping[re.Pattern, Mapping[str, str]] = {
        **generate_ncbi_genbank_nt_regex(("BA", "DF", "DG", "LD"), {
            "DATABASE": "DDBJ",
            "TYPE": "CON division"
        }),
        **generate_ncbi_genbank_nt_regex(("AN",), {
            "DATABASE": "EMBL",
            "TYPE": "CON division"
        }),
        **generate_ncbi_genbank_nt_regex((
            "CH", "CM", "DS", "EM", "EN", "EP", "EQ", "FA", "GG", "GL", "JH", "KB", "KD", "KE", "KI", "KK",
            "KL",
            "KN", "KQ", "KV", "KZ", "ML", "MU"
        ), {
            "DATABASE": "NCBI",
            "TYPE": "CON division"
        }),
        **generate_ncbi_genbank_nt_regex((
            "C", "AT", "AU", "AV", "BB", "BJ", "BP", "BW", "BY", "CI", "CJ", "DA", "DB", "DC", "DK", "FS", "FY",
            "HX", "HY", "LU", "OH"
        ), {
            "DATABASE": "DDBJ",
            "TYPE": "EST"
        }),
        **generate_ncbi_genbank_nt_regex(("F",), {
            "DATABASE": "EMBL",
            "TYPE": "EST"
        }),
        **generate_ncbi_genbank_nt_regex((
            "H", "N", "T", "R", "W", "AA", "AI", "AW", "BE", "BF", "BG", "BI", "BM", "BQ", "BU", "CA", "CB",
            "CD",
            "CF", "CK", "CN", "CO", "CV", "CX", "DN", "DR", "DT", "DV", "DW", "DY", "EB", "EC", "EE",
            "EG", "EH", "EL", "ES", "EV", "EW", "EX", "EY", "FC", "FD", "FE", "FF", "FG", "FK", "FL",
            "GD", "GE", "GH", "GO", "GR", "GT", "GW", "HO", "HS", "JG", "JK", "JZ"
        ), {
            "DATABASE": "GenBank",
            "TYPE": "EST"
        }),
        **generate_ncbi_genbank_nt_regex((
            "D", "AB", "LC"
        ), {
            "DATABASE": "DDBJ",
            "TYPE": "Direct submissions"
        }),
        **generate_ncbi_genbank_nt_regex((
            "V", "X", "Y", "Z", "AJ", "AM", "FM", "FN", "HE", "HF", "HG", "FO", "LK", "LL", "LM", "LN", "LO",
            "LR",
            "LS", "LT", "OA", "OB", "OC", "OD", "OE", "OU", "OV", "OW", "OX", "OY", "OZ",
        ), {
            "DATABASE": "EMBL",
            "TYPE": "Direct submissions"
        }),
        **generate_ncbi_genbank_nt_regex((
            "U", "AF", "AY", "DQ", "EF", "EU", "FJ", "GQ", "GU", "HM", "HQ", "JF", "JN", "JQ", "JX", "KC", "KF",
            "KJ", "KM", "KP", "KR", "KT", "KU", "KX", "KY", "MF", "MG", "MH", "MK", "MN", "MT", "MW", "MZ",
            "OK", "OL", "OM", "ON", "OP", "OQ"
        ), {
            "DATABASE": "GenBank",
            "TYPE": "Direct submissions"
        }),
    }


# AP","BS                                                      DDBJ           Genome project data
# AL","BX","CR","CT","CU","FP","FQ","FR                                    EMBL           Genome project data
# AE","CP","CY                                                   GenBank        Genome project data
# AG","DE","DH","FT","GA","LB                                          DDBJ           GSS
# B","AQ","AZ","BH","BZ","CC","CE","CG","CL","CW","CZ","DU","DX","ED","EI","EJ","EK","         GenBank        GSS
#    ER","ET","FH","FI","GS","HN","HR","JJ","JM","JS","JY","KG","KO","KS","MJ
# AK                                                         DDBJ           cDNA projects
# AC","DP                                                      GenBank        HTGS
# E","BD","DD","DI","DJ","DL","DM","FU","FV","FW","FZ","GB","HV","HW","HZ","LF","LG","         DDBJ           Patents
#    LV","LX","LY","LZ","MA","MB","MC","MD","ME","OF","OG","OI","OJ","PA","PB","PC","
#    PD","PE
# A","AX","CQ","CS","FB","GM","GN","HA","HB","HC","HD","HH","HI","JA","JB","JC","JD","         EMBL           Patents (nucleotide only
#    JE","LP","LQ","MP","MQ","MR","MS
# I","AR","DZ","EA","GC","GP","GV","GX","GY","GZ","HJ","HK","HL","KH","MI","MM","MO","         GenBank        Patents (nucleotide)
#    MV","MX","MY","OO
# G","BV","GF                                                    GenBank        STS
# BR                                                         DDBJ           TPA
# BN                                                         EMBL           TPA
# BK                                                         GenBank        TPA
# HT","HU                                                      DDBJ           TPA CON division
# BL","GJ","GK                                                   GenBank        TPA CON division
# EZ","HP","JI","JL","JO","JP","JR","JT","JU","JV","JW","KA                        GenBank        TSA
# FX","LA","LE","LH","LI","LJ                                          DDBJ           TSA
# S                                                          GenBank        From journal scanning
# AD                                                         GenBank        From GSDB
# AH                                                         GenBank        Segmented set header
# AS                                                         GenBank        Other - not currently being used
# BC                                                         GenBank        MGC project
# BT                                                         GenBank        FLI-cDNA projects
# J","K","L","M                                                    GenBank        From GSDB direct submissions
# N                                                          GenBank/DDBJ   N0-N2 were used intially by both groups but have
#                                                                             been removed from circulation; N2-N9 are ESTs

class AnalysisSetChromosomePlacedGenbankMatcher(AccessionMatcherRuleType):
    _regex = re.compile(r"^chr([0-9]+|X|Y|M)_(.+)v(\d+)?(_alt|_random|_fix)?$")

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
                "TYPE": {
                    "_alt": "Alternate Loci",
                    "_random": "Unlocalized Scaffolds",
                    "_fix": "Fix Patches"
                }.get(groups[3], "N/A")
            }
        )


class AnalysisSetChromosomeUnplacedGenbankMatcher(AccessionMatcherRuleType):
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
                "IS_DECOY": str(groups[2] is not None)
            }
        )


class AnalysisSetChromosomeMatcher(AccessionMatcherRuleType):
    _regex = re.compile(r"^chr([0-9]+|X|Y|M|EBV)$")

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        match_result = self._regex.match(accession)
        if match_result is None:
            return None
        groups = match_result.groups()
        return AccessionMatchResult(
            toplevel="Analysis Set Chromosome",
            details={
                "NUMBER": groups[0],
                "TYPE": {
                    "X": "Sex",
                    "Y": "Sex",
                    "M": "Mitochondria",
                    "EBV": "Epstein-Barr Virus"
                }.get(groups[0], "NORMAL")
            }
        )


class AnalysisSetContigMatcher(ChainAccessionMatcherRuleType):
    _rule_chain: Final[List[AccessionMatcherRuleType]] = [
        AnalysisSetChromosomeMatcher,
        AnalysisSetChromosomeUnplacedGenbankMatcher,
        AnalysisSetChromosomePlacedGenbankMatcher,
        AnalysisSetChromosomeHLA
    ]


class MasterAccessionMatcher(ChainAccessionMatcherRuleType):
    _rule_chain: Final[List[AccessionMatcherRuleType]] = [
        MasterEnsembleIDMatcher,
        NCBIRefSeqMatcher,
        AnalysisSetContigMatcher,
        NCBIGenBankAccessionMatcher
    ]


def infer_accession_type(accession: str) -> Optional[AccessionMatchResult]:
    mam = MasterAccessionMatcher()
    return mam.match(accession)


if __name__ == "__main__":
    print(infer_accession_type("ENSMUSG00000017167.6"))
    print(infer_accession_type("NM_AAAAAAAA.6"))
