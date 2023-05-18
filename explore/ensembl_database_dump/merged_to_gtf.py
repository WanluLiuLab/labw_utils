import pandas as pd

from labw_utils.bioutils.datastructure.gv.gene import DumbGene
from labw_utils.bioutils.parser.gtf import GtfIteratorWriter
from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.record.feature import FeatureInterface, Feature
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List

_lh = get_logger(__name__)
setup_basic_logger()

if __name__ == '__main__':
    _lh.info("Reading data...")
    transcript_table = pd.read_parquet("ensdb_transcripts.parquet.d").dropna(subset=[
        "seqname", "start", "end"
    ])
    exons_table = pd.read_parquet("ensdb_exons.parquet.d").dropna(subset=[
        "seqname", "start", "end"
    ])
    genes_table = pd.read_parquet("ensdb_genes.parquet.d").dropna(subset=[
        "seqname", "start", "end"
    ])

    _lh.info("Merging data...")
    exons_table = exons_table.join(
        transcript_table[
            ["ensdb_pk_transcript_id", "source", "ensdb_transcript_id", "ensdb_transcript_version"]
        ].set_index('ensdb_pk_transcript_id'),
        on="ensdb_pk_transcript_id",
        how="inner"
    )
    transcript_table = transcript_table.join(
        genes_table.drop(columns=['start', 'end', 'strand', 'seqname']).set_index('ensdb_pk_gene_id'),
        on="ensdb_pk_gene_id",
        how="inner"
    )

    gi: List[FeatureInterface] = []
    for exon in tqdm(list(exons_table.itertuples()), desc="Converting exons..."):
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
                "transcript_id": exon.ensdb_transcript_id + "." + str(exon.ensdb_transcript_version),
                "exon_id": exon.ensdb_exon_id + "." + str(exon.ensdb_exon_version),
                "ensdb_exon_created_date": str(exon.ensdb_exon_created_date),
                "ensdb_exon_modified_date": str(exon.ensdb_exon_modified_date),
                "ensdb_pk_exon_id": exon.ensdb_pk_exon_id,
            }
        ))
    for transcript in tqdm(list(transcript_table.itertuples()), desc="Converting transcripts..."):
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
                "transcript_id": transcript.ensdb_transcript_id + "." + str(transcript.ensdb_transcript_version),
                "ensdb_transcript_created_date": str(transcript.ensdb_transcript_created_date),
                "ensdb_transcript_modified_date": str(transcript.ensdb_transcript_modified_date),
                "ensdb_pk_transcript_id": transcript.ensdb_pk_transcript_id,
                "gene_id": transcript.ensdb_gene_id + "." + str(transcript.ensdb_gene_version),
                "ensdb_gene_created_date": str(transcript.ensdb_gene_created_date),
                "ensdb_gene_modified_date": str(transcript.ensdb_gene_modified_date),
                "ensdb_pk_gene_id": transcript.ensdb_pk_gene_id
            }
        ))
    GtfIteratorWriter.write_iterator(
        gi,
        "ensdb_unparsed.gtf"
    )
    GtfIteratorWriter.write_iterator(
        DiploidGeneTree.from_feature_iterator(
            gi,
            gene_implementation=DumbGene
        ).to_feature_iterator(),
        "ensdb.gtf"
    )
