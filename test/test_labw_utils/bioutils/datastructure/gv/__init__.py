from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.record.gtf import parse_record

exon_kwargs = {
    "is_checked": True,
    "shortcut": False
}

exons_str = [
    'chr1\tNA\texon\t5\t10\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.1"; exon_number 1',
    'chr1\tNA\texon\t15\t20\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.2"; exon_number 2',
    'chr1\tNA\texon\t25\t30\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.3"; exon_number 3'
]

exons = list(map(lambda exon_str: Exon(data=parse_record(exon_str), **exon_kwargs), exons_str))
