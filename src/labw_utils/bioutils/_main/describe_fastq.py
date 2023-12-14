"""
describe_fastq.py -- Lite Python-implemented fastqc.

.. versionadded:: 1.0.2
.. versionchanged:: 1.0.3
    The param --out is truly optional.
"""
import argparse
import os

from labw_utils import UnmetDependenciesError
from labw_utils.bioutils.parser.fasta import FastaIterator
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.shutil_helper import wc_c
from labw_utils.typing_importer import List, Optional, Literal

try:
    import numpy as np
except ImportError:
    raise UnmetDependenciesError("numpy")

from labw_utils.bioutils.algorithm.sequence import get_gc_percent, decode_phred33
from labw_utils.bioutils.parser.fastq import FastqIterator
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.appender import load_table_appender_class, TableAppenderConfig


_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m labw_utils.bioutils transcribe",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path of Input FASTA & FASTQ",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=False,
        help="Path of output directory",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--with_per_pos_base_qual",
        help="Store extension_stats.tsv. Recommended for NGS but not for TGS.",
        action="store_true",
    )
    parser.add_argument(
        "--full_header",
        help="Use fill header.",
        action="store_true",
    )
    parser.add_argument(
        "--input_fmt",
        help="Format of input",
        choices=["fastq", "fasta"],
        required=True,
    )
    return parser


def fastqc(
    filepath: str,
    outdir_path: Optional[str] = None,
    with_per_pos_base_qual: bool = True,
    input_fmt: Literal["fastq", "fasta"] = "fastq",
    full_header: bool = True,
):
    _lh.info("Start parsing '%s'...", filepath)

    if not os.path.exists(filepath):
        _lh.error("File '%s' not exist!", filepath)
    outdir_path = filepath + ".stats.d" if outdir_path is None else outdir_path
    os.makedirs(outdir_path, exist_ok=True)

    if input_fmt == "fasta" and with_per_pos_base_qual:
        with_per_pos_base_qual = False
        _lh.warning("with_per_pos_base_qual ias set to FALSE since it is not possible in FASTA input")

    max_read_length = 0
    with load_table_appender_class("TSVTableAppender")(
        os.path.join(outdir_path, "all"),
        ("SEQID", "GC", "LEN", "MEANQUAL"),
        TableAppenderConfig(),
    ) as appender:
        if input_fmt == "fastq":
            for fastq_record in FastqIterator(
                filename=filepath,
                show_tqdm=(wc_c(filepath) < 10 * (1 << 20)),
                full_header=full_header,
            ):
                len_record = len(fastq_record)
                appender.append(
                    (
                        fastq_record.seq_id,
                        get_gc_percent(fastq_record.sequence),
                        len_record,
                        np.mean(list(decode_phred33(fastq_record.quality))),
                    )
                )
                max_read_length = max(max_read_length, len_record)
        else:
            for fasta_record in FastaIterator(
                filename=filepath,
                show_tqdm=(wc_c(filepath) < (1 << 20)),
                full_header=full_header,
            ):
                len_record = len(fasta_record)
                appender.append(
                    (
                        fasta_record.seq_id,
                        get_gc_percent(fasta_record.sequence),
                        len_record,
                        0,
                    )
                )
    if with_per_pos_base_qual:
        quality_sum = np.zeros(max_read_length, dtype=np.double)
        quality_cnt = np.zeros(max_read_length, dtype=np.double)
        for fastq_record in FastqIterator(filename=filepath):
            current_quality = np.array(list(decode_phred33(fastq_record.quality)))
            quality_sum[0 : current_quality.shape[0]] += current_quality
            quality_cnt[0 : current_quality.shape[0]] += 1
        quality = quality_sum / quality_cnt
        with load_table_appender_class("TSVTableAppender")(
            os.path.join(outdir_path, "extension_stat"),
            ("POS", "QUAL"),
            TableAppenderConfig(),
        ) as appender:
            for pos, qual in enumerate(quality):
                appender.append((pos, qual))


def main(args: List[str]):
    argv = create_parser().parse_args(args)
    fastqc(argv.input, argv.outdir, argv.with_per_pos_base_qual, argv.input_fmt, argv.full_header)
