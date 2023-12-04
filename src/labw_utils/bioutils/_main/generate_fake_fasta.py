"""
generate_fake_fasta.py -- Create fake organism.

.. versionadded:: 1.0.2
"""
__all__ = (
    "create_parser",
    "main",
)

import argparse
import random

from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.typing_importer import List


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m labw_utils.bioutils generate_fake_fasta",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Path to output FASTA",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--chr_len",
        required=False,
        help="Length of each chromosome",
        nargs="?",
        type=int,
        action="store",
        default=10000,
    )
    parser.add_argument(
        "--chr_num",
        required=False,
        help="Number of chromosome",
        nargs="?",
        type=int,
        action="store",
        default=5,
    )
    parser.add_argument(
        "--n_3p",
        required=False,
        help="Proportion of 3' Ns",
        nargs="?",
        type=float,
        action="store",
        default=0.1,
    )
    parser.add_argument(
        "--n_5p",
        required=False,
        help="Proportion of 5' Ns",
        nargs="?",
        type=float,
        action="store",
        default=0.1,
    )
    parser.add_argument(
        "--gc",
        required=False,
        help="GC Percent",
        nargs="?",
        type=float,
        action="store",
        default=0.4,
    )
    parser.add_argument(
        "--centromere_offset",
        required=False,
        help="Offset of centromere",
        nargs="?",
        type=float,
        action="store",
        default=0.3,
    )
    parser.add_argument(
        "--centromere_len",
        required=False,
        help="Length of centromere",
        nargs="?",
        type=float,
        action="store",
        default=0.2,
    )
    parser.add_argument(
        "--split_at",
        required=False,
        help="Length of split, 0 for unsplit",
        nargs="?",
        type=int,
        action="store",
        default=80,
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    with get_writer(args.out) as writer:
        for seqid in range(args.chr_num):
            chr_contents = list(
                map(
                    lambda _: random.choice(random.choices(("CG", "AT"), (args.gc, 1 - args.gc))[0]),
                    range(args.chr_len),
                )
            )

            def mask(start: float, end: float):
                start, end = int(start), int(end)
                chr_contents[start:end] = ["N"] * int(end - start)

            if args.n_5p != 0:
                mask(0, args.chr_len * args.n_5p)
            if args.n_3p != 0:
                mask(args.chr_len - args.chr_len * args.n_3p, args.chr_len)
            if args.centromere_len != 0:
                mask(
                    args.chr_len * args.centromere_offset,
                    args.chr_len * (args.centromere_offset + args.centromere_len),
                )
            chr_contents = "".join(chr_contents)
            if args.split_at != 0:
                chr_contents = "\n".join(
                    list(chr_contents[i : i + args.split_at] for i in range(0, len(chr_contents), args.split_at))
                )
            writer.write(f">{seqid}\n{chr_contents}\n")
