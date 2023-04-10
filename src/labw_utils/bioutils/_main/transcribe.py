"""
transcribe.py -- General-purposed stranded transcription, from reference genome to reference cDNA.
"""

__all__ = (
    "create_parser",
    "main"
)

import argparse
from typing import List

from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gv_helper import transcribe


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m labw_utils.bioutils transcribe",
        description=__doc__.splitlines()[1]
    )
    parser = FrontendOptSpecs.patch(parser, "-f")
    parser = FrontendOptSpecs.patch(parser, "-g")
    parser.add_argument(
        '-o', '--out',
        required=True,
        help="Path of Output cDNA FASTA",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--no_write_single_transcript',
        help="Stop splitting cDNA of each isoform into separate file",
        action='store_true'
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    fv = FastaViewFactory(args.fasta)
    transcribe(
        gv,
        args.out,
        fv,
        write_single_transcript=not args.no_write_single_transcript
    )
