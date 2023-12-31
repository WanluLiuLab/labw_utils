import os
import tempfile

from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.commonutils.lwio.safe_io import get_writer

gene_gtf = """
chrI	ncbiRefSeq	exon	4221	4358	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "1"; exon_id "NM_058260.4.1"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	5195	5296	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "2"; exon_id "NM_058260.4.2"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	6037	6327	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "3"; exon_id "NM_058260.4.3"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	9727	9846	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "4"; exon_id "NM_058260.4.4"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	10095	10148	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "5"; exon_id "NM_058260.4.5"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	11641	11689	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "1"; exon_id "NM_058259.4.1"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	11641	11689	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "1"; exon_id "NM_058259.4.1"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	16473	16585	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "3"; exon_id "NM_058259.4.3"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	14951	15160	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "2"; exon_id "NM_058259.4.2"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	15103	15160	.	+	.	gene_id "nlp-40"; transcript_id "NM_001306277.1"; exon_number "1"; exon_id "NM_001306277.1.1"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	16473	16585	.	+	.	gene_id "nlp-40"; transcript_id "NM_001306277.1"; exon_number "2"; exon_id "NM_001306277.1.2"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	8484719	8484913	.	-	.	gene_id "D1081.6"; transcript_id "NM_059899.3"; exon_number "1"; exon_id "NM_059899.3.1"; gene_name "D1081.6";
chrI	ncbiRefSeq	transcript	6492885	6493553	0	-	.	gene_id "mdt-18"; transcript_id "NM_001322685.1"; gene_name "mdt-18"; 
"""


def test_gene() -> None:
    with tempfile.TemporaryDirectory() as test_path:
        with get_writer(os.path.join(test_path, "1.gtf")) as writer:
            writer.write(gene_gtf)
        gv = GeneViewFactory.from_file(os.path.join(test_path, "1.gtf"))
        assert list(gv.iter_gene_ids()) == ["homt-1", "nlp-40", "D1081.6", "mdt-18"]
        assert list(gv.iter_transcript_ids()) == [
            "NM_058260.4",
            "NM_058259.4",
            "NM_001306277.1",
            "NM_059899.3",
            "NM_001322685.1",
        ]
        assert gv.get_transcript("NM_058260.4").get_nth_exon(0).start == 4221
        gv.standardize()
        assert list(gv.iter_transcript_ids()) == ["NM_058260.4", "NM_058259.4", "NM_001306277.1", "NM_059899.3"]
        gv.to_file(os.path.join(test_path, "3.gtf"))
        os.system(f"gedit {test_path}/*.gtf")
