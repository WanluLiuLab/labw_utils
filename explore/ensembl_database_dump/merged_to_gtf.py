import pandas as pd

from labw_utils.bioutils.parser.gtf import GtfIteratorWriter
from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.record.feature import FeatureInterface, Feature
from labw_utils.typing_importer import List

if __name__ == '__main__':

    transcript_table = pd.read_parquet("ensdb_transcripts.parquet.d").dropna(subset=[
        "seqname", "start", "end"
    ])
    exons_table = pd.read_parquet("ensdb_exons.parquet.d").dropna(subset=[
        "seqname", "start", "end"
    ])
    genes_table = pd.read_parquet("ensdb_genes.parquet.d").dropna(subset=[
        "seqname", "start", "end"
    ])
    gi: List[FeatureInterface] = []
    exons_table = exons_table.join(
        transcript_table[
            ["ensdb_pk_transcript_id", "source", "ensdb_transcript_id", "ensdb_transcript_version"]
        ].set_index('ensdb_pk_transcript_id'),
        on="ensdb_pk_transcript_id",
        how="inner"
    )

    for exon in exons_table.itertuples():
        if exon.strand == 1:
            strand = "+"
        elif exon.strand == -1:
            strand = "-"
        else:
            strand = "."
        gi.append(Feature(
            seqname=exon.seqname,
            source=exon.source,
            start=exon.start,
            end=exon.end,
            strand=strand,
            feature="exon",
            score=None,
            frame=None,
            attribute={
                "transcript_id": exon.ensdb_transcript_id+"."+str(exon.ensdb_transcript_version),
                "exon_id":exon.ensdb_exon_id+"."+str(exon.ensdb_exon_version),
                "ensdb_exon_created_date": exon.ensdb_exon_created_date,
                "ensdb_exon_modified_date": exon.ensdb_exon_modified_date,
                "ensdb_pk_exon_id": exon.ensdb_pk_exon_id,
            }
        ))
        print(gi[-1])
    for transcript in transcript_table.itertuples():
        if transcript.strand == 1:
            strand = "+"
        elif transcript.strand == -1:
            strand = "-"
        else:
            strand = "."
        gi.append(Feature(
            seqname=transcript.seqname,
            source=transcript.source,
            start=transcript.start,
            end=transcript.end,
            strand=strand,
            feature="transcript",
            score=None,
            frame=None,
            attribute={
                "transcript_id": exon.ensdb_transcript_id+"."+str(exon.ensdb_transcript_version),
                "ensdb_transcript_created_date": exon.ensdb_transcript_created_date,
                "ensdb_transcript_modified_date": exon.ensdb_transcript_modified_date,
                "ensdb_pk_transcript_id": exon.ensdb_pk_transcript_id
            }
        ))
        print(gi[-1])

    GtfIteratorWriter.write_iterator(
        DiploidGeneTree.from_feature_iterator(gi).to_feature_iterator(),
        "ensdb.gtf"
    )
