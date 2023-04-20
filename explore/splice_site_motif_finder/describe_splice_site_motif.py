import pandas as pd

from labw_utils.bioutils.algorithm.sequence import reverse_complement
from labw_utils.bioutils.datastructure import gene_tree
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.importer.tqdm_importer import tqdm

if __name__ == "__main__":
    gt= gene_tree.GeneTree.from_gtf_file("ce11.ncbiRefSeq.gtf")
    fav = FastaViewFactory("ce11.fa")
    starts = pd.DataFrame(
        data=0,
        columns=["A", "G", "C", "T", "N"],
        index=range(-100, 100)
    )
    ends = pd.DataFrame(
        data=0,
        columns=["A", "G", "C", "T", "N"],
        index=range(-100, 100)
    )
    for transcript in tqdm(gt.transcript_values, total=gt.number_of_transcripts):
        seqname = transcript.seqname
        strand = transcript.strand
        for splice_site_start, splice_site_end in transcript.splice_sites:
            start_seq = fav.sequence(seqname, splice_site_start - 100, splice_site_start + 100)
            end_seq = fav.sequence(seqname, splice_site_end - 100, splice_site_end + 100)

            if not transcript.strand:
                start_seq = reverse_complement(start_seq)
                end_seq = reverse_complement(end_seq)
            for i, pos in enumerate(start_seq, start = -100):
                pos = pos.upper()
                if pos not in "AGCT":
                    pos = "N"
                starts.loc[i, pos] += 1
            for i, pos in enumerate(end_seq, start = -100):
                pos = pos.upper()
                if pos not in "AGCT":
                    pos = "N"
                ends.loc[i, pos] += 1


    starts.to_csv("start.csv")
    ends.to_csv("ends.csv")
