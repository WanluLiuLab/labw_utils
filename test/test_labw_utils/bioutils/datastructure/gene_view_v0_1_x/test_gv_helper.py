import os
import tempfile

import pytest

import labw_utils.bioutils.datastructure.gene_view_v0_1_x.gv_helper as gvh
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.commonutils.lwio.safe_io import get_writer

gene_gtf = """
1	StringTie	transcript	14361	29359	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; cov "14.715067"; FPKM "1.148767"; TPM "1.453300";
1	StringTie	exon	14361	14829	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; cov "13.400369";
1	StringTie	exon	14970	15038	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "2"; cov "7.045532";
1	StringTie	exon	15796	16765	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "3"; cov "3.927150";
1	StringTie	exon	16858	17055	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "4"; cov "27.734772";
1	StringTie	exon	17233	17368	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "5"; cov "26.952866";
1	StringTie	exon	17606	17742	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "6"; cov "29.338306";
1	StringTie	exon	17915	18061	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "7"; cov "26.099281";
1	StringTie	exon	18268	18366	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "8"; cov "26.265776";
1	StringTie	exon	24738	24891	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "9"; cov "31.549301";
1	StringTie	exon	29321	29359	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "10"; cov "13.560440";
1	StringTie	transcript	14361	25012	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; cov "1.460519"; FPKM "0.114019"; TPM "0.144245";
1	StringTie	exon	14361	14829	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "1"; cov "0.432010";
1	StringTie	exon	14970	15038	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "2"; cov "1.040070";
1	StringTie	exon	15796	16765	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "3"; cov "1.043148";
1	StringTie	exon	16858	17055	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "4"; cov "2.787202";
1	StringTie	exon	17233	17368	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "5"; cov "2.707201";
1	StringTie	exon	17606	17742	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "6"; cov "2.793163";
1	StringTie	exon	17915	18061	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "7"; cov "2.173379";
1	StringTie	exon	24738	25012	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "8"; cov "2.175562";
"""

gene_gtf_with_duplicated_transcripts_in_one_gene = """
1	StringTie	transcript	14361	29359	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; cov "14.715067"; FPKM "1.148767"; TPM "1.453300";
1	StringTie	exon	14361	14829	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; cov "13.400369";
1	StringTie	exon	14970	15038	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "2"; cov "7.045532";
1	StringTie	exon	15796	16765	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "3"; cov "3.927150";
1	StringTie	exon	16858	17055	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "4"; cov "27.734772";
1	StringTie	exon	17233	17368	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "5"; cov "26.952866";
1	StringTie	exon	17606	17742	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "6"; cov "29.338306";
1	StringTie	exon	17915	18061	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "7"; cov "26.099281";
1	StringTie	exon	18268	18366	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "8"; cov "26.265776";
1	StringTie	exon	24738	24891	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "9"; cov "31.549301";
1	StringTie	exon	29321	29359	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "10"; cov "13.560440";

1	StringTie	transcript	14361	25012	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; cov "1.460519"; FPKM "0.114019"; TPM "0.144245";
1	StringTie	exon	14361	14829	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "1"; cov "0.432010";
1	StringTie	exon	14970	15038	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "2"; cov "1.040070";
1	StringTie	exon	15796	16765	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "3"; cov "1.043148";
1	StringTie	exon	16858	17055	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "4"; cov "2.787202";
1	StringTie	exon	17233	17368	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "5"; cov "2.707201";
1	StringTie	exon	17606	17742	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "6"; cov "2.793163";
1	StringTie	exon	17915	18061	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "7"; cov "2.173379";
1	StringTie	exon	24738	25012	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "8"; cov "2.175562";

1	StringTie	transcript	14361	29359	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; cov "14.715067"; FPKM "1.148767"; TPM "1.453300";
1	StringTie	exon	14361	14829	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "1"; cov "13.400369";
1	StringTie	exon	14970	15038	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "2"; cov "7.045532";
1	StringTie	exon	15796	16765	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "3"; cov "3.927150";
1	StringTie	exon	16858	17055	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "4"; cov "27.734772";
1	StringTie	exon	17233	17368	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "5"; cov "26.952866";
1	StringTie	exon	17606	17742	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "6"; cov "29.338306";
1	StringTie	exon	17915	18061	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "7"; cov "26.099281";
1	StringTie	exon	18268	18366	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "8"; cov "26.265776";
1	StringTie	exon	24738	24891	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "9"; cov "31.549301";
1	StringTie	exon	29321	29359	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.3"; exon_number "10"; cov "13.560440";

1	StringTie	transcript	14361	25012	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; cov "1.460519"; FPKM "0.114019"; TPM "0.144245";
1	StringTie	exon	14361	14829	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; exon_number "1"; cov "0.432010";
1	StringTie	exon	14970	15038	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; exon_number "2"; cov "1.040070";
1	StringTie	exon	15796	16765	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; exon_number "3"; cov "1.043148";
1	StringTie	exon	16858	17055	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; exon_number "4"; cov "2.787202";
1	StringTie	exon	17233	17368	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; exon_number "5"; cov "2.707201";
1	StringTie	exon	17606	17742	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; exon_number "6"; cov "2.793163";
1	StringTie	exon	17915	18061	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; exon_number "7"; cov "2.173379";
1	StringTie	exon	24738	25010	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.4"; exon_number "8"; cov "2.175562";
"""

gene_gtf_for_exon_supersets = """
1	StringTie	transcript	14361	29359	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; cov "14.715067"; FPKM "1.148767"; TPM "1.453300";
1	StringTie	exon	14361	14829	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; cov "13.400369";
1	StringTie	exon	14970	15038	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "2"; cov "7.045532";
1	StringTie	exon	15796	16765	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "3"; cov "3.927150";
1	StringTie	exon	16858	17055	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "4"; cov "27.734772";
1	StringTie	exon	17233	17368	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "5"; cov "26.952866";
1	StringTie	exon	17606	17742	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "6"; cov "29.338306";
1	StringTie	exon	17915	18061	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "7"; cov "26.099281";
1	StringTie	exon	18268	18366	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "8"; cov "26.265776";
1	StringTie	exon	24738	24891	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "9"; cov "31.549301";
1	StringTie	exon	29321	29359	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "10"; cov "13.560440";
1	StringTie	transcript	14361	25012	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; cov "1.460519"; FPKM "0.114019"; TPM "0.144245";
1	StringTie	exon	14361	14829	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "1"; cov "0.432010";
1	StringTie	exon	14970	15038	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "2"; cov "1.040070";
1	StringTie	exon	15796	16765	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "3"; cov "1.043148";
1	StringTie	exon	16858	17055	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "4"; cov "2.787202";
1	StringTie	exon	17233	17368	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "5"; cov "2.707201";
1	StringTie	exon	17606	17742	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "6"; cov "2.793163";
1	StringTie	exon	17915	18061	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "7"; cov "2.173379";
1	StringTie	exon	24738	25012	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "8"; cov "2.175562";
"""


@pytest.fixture(scope="module", autouse=False)
def initialize_module() -> str:
    """
    This function sets up a directory for testing
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        with get_writer(os.path.join(tmpdir, "1.gtf.gz")) as fh:
            fh.write(gene_gtf)
        with get_writer(os.path.join(tmpdir, "gene_gtf_with_duplicated_transcripts_in_one_gene.gtf")) as fh:
            fh.write(gene_gtf_with_duplicated_transcripts_in_one_gene)
        with get_writer(os.path.join(tmpdir, "gene_gtf_for_exon_supersets.gtf")) as fh:
            fh.write(gene_gtf_for_exon_supersets)
        yield tmpdir


def test_gene_gtf_with_duplicated_transcripts_in_one_gene(initialize_module):
    file_name = os.path.join(initialize_module, "gene_gtf_with_duplicated_transcripts_in_one_gene.gtf")
    gv = GeneViewFactory.from_file(file_name)
    gvh.gv_dedup(gv, by_splice_site=True, assume_no_cross_gene_duplication=True)
    assert gv.number_of_transcripts == 2


def test_gene_gtf_with_duplicated_transcripts_in_one_gene_by_splice_sites(initialize_module):
    file_name = os.path.join(initialize_module, "gene_gtf_with_duplicated_transcripts_in_one_gene.gtf")
    gv = GeneViewFactory.from_file(file_name)
    gvh.gv_dedup(gv, by_splice_site=False, assume_no_cross_gene_duplication=True)
    assert gv.number_of_transcripts == 3
