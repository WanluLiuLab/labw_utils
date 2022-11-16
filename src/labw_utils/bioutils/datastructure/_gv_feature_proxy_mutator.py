import bisect
import math

from labw_utils.bioutils.datastructure._gv_errors import *
from labw_utils.bioutils.datastructure.gv_feature_proxy import DEFAULT_SORT_EXON_EXON_STRAND_POLICY, Transcript, Exon, \
    Gene



class GeneMutator:

    @staticmethod
    def del_transcript(gene: Gene, transcript_id: str):
        index = gene._transcript_ids.index(transcript_id)
        gene._transcript_ids.pop(index)
        gene._transcripts.pop(index)

    @staticmethod
    def replace_transcript(
            gene: Gene,
            transcript: Transcript,
            **kwargs
    ):
        gene.del_transcript(transcript.transcript_id)
        gene.add_transcript(transcript, **kwargs)

    @staticmethod
    def fast_add_transcript(
            gene: Gene,
            transcript: Transcript,
    ):
        GeneMutator.add_transcript(gene, transcript, False, False)

    @staticmethod
    def add_transcript(
            gene: Gene,
            transcript: Transcript,
            check_same_chrome: bool = True,
            check_same_strand: bool = True
    ):
        if transcript.transcript_id in gene._transcript_ids:
            raise DuplicatedTranscriptIDError(
                f"Transcript ID {transcript.transcript_id} duplicated"
            )
        new_pos = bisect.bisect_left(gene._transcripts, transcript)
        if check_same_chrome and transcript.seqname != gene.seqname:
            raise TranscriptInAGeneOnDifferentChromosomeError(
                f"gene.seqname={gene.seqname}, while transcript.seqname={transcript.seqname}"
            )
        if check_same_strand and transcript.strand != gene.strand and transcript.strand != "." and gene.strand != ".":
            raise TranscriptInAGeneOnDifferentStrandError
        gene._transcript_ids.insert(new_pos, transcript.transcript_id)
        gene._transcripts.insert(new_pos, transcript)
