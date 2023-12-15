"""
transcribe.py -- General-purposed stranded transcription, from reference genome to reference cDNA.

.. versionadded:: 1.0.2
"""

__all__ = (
    "create_parser",
    "main",
)

import argparse

from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.datastructure.gene_tree_helper import transcribe
from labw_utils.bioutils.datastructure.gv.gene import DumbGene

from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.typing_importer import List


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m labw_utils.bioutils transcribe",
        description=__doc__.splitlines()[1],
    )
    parser = FrontendOptSpecs.patch(parser, "-f")
    parser = FrontendOptSpecs.patch(parser, "-g")
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Path of Output cDNA FASTA",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--no_write_single_transcript",
        help="Stop splitting cDNA of each isoform into separate file",
        action="store_true",
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    gv = DiploidGeneTree.from_gtf_file(args.gtf, gene_implementation=DumbGene)
    fv = FastaViewFactory(args.fasta)
    transcribe(
        gv,
        args.out,
        fv,
        write_single_transcript=not args.no_write_single_transcript,
    )
