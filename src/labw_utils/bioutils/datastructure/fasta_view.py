"""
``labw_utils.bioutils.datastructure.fasta_view`` -- General FASTA reader

Can provide random access to FASTA files, compressed or non-compressed.

Highlights: This utility can read all format supported by :py:mod:`commonutils.lwio`,
while others require Block GZipped ones.

.. note::
    Although this module supports all format supported by :py:mod:`commonutils.lwio`,
    it is recommended for user to compress their files using ``bgzip`` and index them using ``tabix``.

.. warning::
    This module uses 0-based ``[)`` indexing!

.. versionadded:: 1.0.2
"""

__all__ = (
    'FastaViewType',
    'FastaViewFactory',
    "SeekTooFarError",
    "ChromosomeNotFoundError",
    "FromGreaterThanToError",
    "FastaViewInvalidRegionError",
    "DuplicatedChromosomeNameError",
    "split_fasta"
)

import functools
import os
from abc import abstractmethod, ABC

from labw_utils.bioutils.datastructure.fai_view import FastaIndexView
from labw_utils.bioutils.parser.fai import FastaIndexNotWritableError
from labw_utils.bioutils.parser.fasta import FastaIterator
from labw_utils.commonutils.lwio.file_system import file_exists
from labw_utils.commonutils.lwio.safe_io import get_reader, get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.shutil_helper import wc_c
from labw_utils.typing_importer import List, Union, Tuple, Dict, Optional, IO, Iterable

_lh = get_logger(__name__)

QueryTupleType = Union[Tuple[str, int, int], Tuple[str, int], Tuple[str]]
FASTA_SPLIT_SEQNAME_OPTIONS = (
    "error",
    "convert",
    "skip"
)


class FastaViewError(ValueError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    pass


class DuplicatedChromosomeNameError(FastaViewError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    def __init__(self, name: str):
        super().__init__(f"Chromosome name {name} duplicated")


class FastaViewInvalidRegionError(FastaViewError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    pass


class SeekTooFarError(FastaViewInvalidRegionError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    def __init__(self, chromosome: str, pos: int, chr_len: int):
        super().__init__(f"Seek {pos}@{chromosome} too far, valid is -1, [0, {chr_len})")


class ChromosomeNotFoundError(FastaViewInvalidRegionError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    def __init__(self, chromosome: str):
        super().__init__(f"Requested chromosome '{chromosome}' not found")


class FromGreaterThanToError(FastaViewInvalidRegionError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    def __init__(self, from_pos: int, to_pos: int):
        super().__init__(f"Requested from_pos {from_pos} > to_pos {to_pos} not allowed!")


class FastaViewType:
    """
    Abstract class of factories.

    .. versionadded:: 1.0.2
    """
    filename: str
    """
    Filename to read
    """

    full_header: bool
    """
    Whether to read in all header.
    If False (default), will read until seeing space or tab.
    See :py:mod:`pybedtools` for more details.
    """

    @abstractmethod
    def sequence(self, chromosome: str, from_pos: int = 0, to_pos: int = -1) -> str:
        """
        Get sequence from FASTA with 0-based [) indexes.

        To read until end, use -1.

        :param chromosome: Chromosome name
        :param from_pos: From which position
        :param to_pos: To which position, use -1 for end
        """
        raise NotImplementedError

    @abstractmethod
    def get_chr_length(self, chromosome: str) -> int:
        """
        Get length of specific chromosome.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def chr_names(self) -> List[str]:
        """
        get chromosome names
        """
        raise NotImplementedError

    def __str__(self) -> str:
        return repr(self)

    @abstractmethod
    def is_valid_region(self, chromosome: str, from_pos: int, to_pos: int):
        """
        Whether a region is valid. See :py:func:`sequence` for details.

        :raises SeekTooFarError: Raise this error if region is not valid.
        :raises ChromosomeNotFoundError: Raise this error if region is not valid.
        """
        raise NotImplementedError

    @abstractmethod
    def close(self):
        """
        Safely close a Fasta.
        """
        raise NotImplementedError

    @abstractmethod
    def __len__(self) -> int:
        raise NotImplementedError

    @abstractmethod
    def __repr__(self):
        raise NotImplementedError

    @abstractmethod
    def __del__(self):
        raise NotImplementedError

    @abstractmethod
    def to_file(self, output_filename: str):
        """
        Write content of this FASTA view to file
        """
        raise NotImplementedError

    @abstractmethod
    def query(self, query: QueryTupleType) -> str:
        """
        :py:func:`sequence` with another interface.
        """
        raise NotImplementedError

    @abstractmethod
    def subset_to_file(
            self,
            output_filename: str,
            querys: Iterable[QueryTupleType],
            output_chr_names: Optional[Iterable[str]] = None
    ):
        raise NotImplementedError

    @abstractmethod
    def __enter__(self):
        raise NotImplementedError

    @abstractmethod
    def __exit__(self, exc_type, exc_val, exc_tb):
        raise NotImplementedError

    @abstractmethod
    def legalize_region_best_effort(self, chromosome: str, from_pos: int = 0, to_pos: int = -1) -> QueryTupleType:
        """
        # TODO
        """
        raise NotImplementedError


class _BaseFastaView(FastaViewType, ABC):
    """
    Base class of other backends.
    """

    def __init__(self, filename: str, full_header: bool = False):
        self.full_header = full_header
        self.filename = filename

    def is_valid_region(self, chromosome: str, from_pos: int, to_pos: int):
        if chromosome not in self.chr_names:
            raise ChromosomeNotFoundError(chromosome)
        chr_len = self.get_chr_length(chromosome)
        if from_pos < 0 or from_pos > chr_len:
            raise SeekTooFarError(chromosome, from_pos, chr_len)
        if to_pos != -1 and to_pos < 0 or to_pos > chr_len:
            raise SeekTooFarError(chromosome, to_pos, chr_len)
        if to_pos != -1 and from_pos > to_pos:
            raise FromGreaterThanToError(from_pos, to_pos)

    def __len__(self) -> int:
        return len(self.chr_names)

    def __repr__(self):
        try:
            return f"Fasta from {self.filename}, len={len(self)}"
        except AttributeError:
            return "Fasta being constructed"

    def __del__(self):
        self.close()

    def to_file(self, output_filename: str):
        with get_writer(output_filename) as writer:
            for k in self.chr_names:
                fa_str = f">{k}\n{self.sequence(k)}\n"
                writer.write(fa_str)

    def query(self, query: QueryTupleType) -> str:
        return self.sequence(*query)

    def subset_to_file(
            self,
            output_filename: str,
            querys: Iterable[QueryTupleType],
            output_chr_names: Optional[Iterable[str]] = None
    ):
        querys = list(querys)
        if output_chr_names is None:
            output_chr_names = list(map(lambda x: x[0], querys))
        else:
            output_chr_names = list(output_chr_names)
        for output_chr_name in output_chr_names:
            if output_chr_names.count(output_chr_name) > 1:
                raise DuplicatedChromosomeNameError(output_chr_name)
        with get_writer(output_filename) as writer:
            for output_chr_name, query in zip(output_chr_names, querys):
                fa_str = f">{output_chr_name}\n{self.query(query)}\n"
                writer.write(fa_str)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def legalize_region_best_effort(self, chromosome: str, from_pos: int = 0, to_pos: int = -1) -> QueryTupleType:
        if chromosome not in self.chr_names:
            raise ChromosomeNotFoundError(chromosome)
        chr_len = self.get_chr_length(chromosome)
        if from_pos < 0:
            from_pos = 0
        elif from_pos > chr_len:
            from_pos = chr_len
        if to_pos == -1:
            to_pos = chr_len
        elif to_pos < 0:
            to_pos = 0
        elif to_pos > chr_len:
            to_pos = chr_len
        return chromosome, from_pos, to_pos


class _MemoryAccessFastaView(_BaseFastaView):
    """
    Fasta whose sequences are read into memory.
    Extremely fast but need lots of memory. Suitable for small files.

    .. versionadded:: 1.0.2
    """

    _all_dict: Dict[str, str]
    """
    Dict[chromosome_name, sequence]
    """

    @property
    def chr_names(self) -> List[str]:
        return list(self._all_dict.keys())

    @functools.lru_cache(maxsize=256, typed=True)
    def get_chr_length(self, chromosome: str) -> int:
        return len(self._all_dict[chromosome])

    def close(self):
        """Nothing to close as is in memory"""
        pass

    def __init__(self, filename: str, full_header: bool = False, show_tqdm: bool = True):
        super().__init__(filename, full_header)
        self._all_dict = {}  # For in-memory reader, will read in all sequences
        self._read_into_mem(show_tqdm=show_tqdm)

    def _read_into_mem(self, show_tqdm: bool) -> None:
        """
        Read FASTA into memory
        """
        for fasta_record in FastaIterator(
                show_tqdm=show_tqdm,
                filename=self.filename,
                full_header=self.full_header
        ):
            if fasta_record.seq_id in self._all_dict:
                raise DuplicatedChromosomeNameError(fasta_record.seq_id)
            self._all_dict[fasta_record.seq_id] = fasta_record.sequence

    def sequence(self, chromosome: str, from_pos: int = 0, to_pos: int = -1):
        self.is_valid_region(chromosome, from_pos, to_pos)
        if to_pos == -1:
            to_pos = self.get_chr_length(chromosome)
        return self._all_dict[chromosome][from_pos:to_pos]


class _DiskAccessFastaView(_BaseFastaView):
    """
    Fasta whose sequence is NOT read into memory.
    Slow but memory-efficient.

    .. versionadded:: 1.0.2
    """

    _fd: IO
    """
    Underlying file descriptor
    """

    _fai: FastaIndexView

    @functools.lru_cache(maxsize=256, typed=True)
    def get_chr_length(self, chromosome: str) -> int:
        return self._fai[chromosome].length

    @property
    def chr_names(self) -> List[str]:
        return list(self._fai.keys())

    def __init__(
            self,
            filename: str,
            full_header: bool = False,
            show_tqdm: bool = True
    ):
        super().__init__(filename, full_header)
        # If has prebuilt index file, read it
        self._fd = get_reader(self.filename)
        index_filename = self.filename + ".fai"
        if file_exists(index_filename):
            if os.path.getmtime(index_filename) - os.path.getmtime(filename) < 0:
                _lh.warning("Index filename %s older than %s, will be rebuilt", index_filename, filename)
            else:
                self._fai = FastaIndexView.from_fai(index_filename, show_tqdm=show_tqdm)
                return
        else:
            _lh.warning("Index filename %s not exist", index_filename)
        self._fai = FastaIndexView.from_fasta(
            filename=filename,
            full_header=full_header,
            show_tqdm=show_tqdm
        )
        try:
            self._fai.write(index_filename)
        except FastaIndexNotWritableError as e:
            _lh.error("Fasta index generated but not writable %s", e)

    def sequence(self, chromosome: str, from_pos: int = 0, to_pos: int = -1) -> str:
        self.is_valid_region(chromosome, from_pos, to_pos)
        chr_fai = self._fai[chromosome]
        """FAI record of this chromosome"""

        if to_pos == -1:
            to_pos = self.get_chr_length(chromosome)

        len_newline = chr_fai.line_len - chr_fai.line_blen
        """Length of newlines"""

        self._fd.seek(chr_fai.offset + from_pos // chr_fai.line_blen * len_newline + from_pos)
        prev_resid = from_pos % chr_fai.line_blen
        """Previous residue. Where the reader is on"""

        lines_to_read = (to_pos - from_pos + prev_resid) // chr_fai.line_blen
        """
        How many full-length line to read
        """

        if (to_pos % chr_fai.line_blen) == 0:
            lines_to_read -= 1

        rets = self._fd.read(
            to_pos - from_pos + lines_to_read * len_newline
        ).replace('\n', '').replace('\r', '')
        return rets

    def close(self):
        try:
            self._fd.close()
        except AttributeError:
            pass


def FastaViewFactory(
        filename: str,
        full_header: bool = False,
        read_into_memory: Optional[bool] = None,
        show_tqdm: bool = True
) -> FastaViewType:
    """
    Initialize a _DiskFasta interface using multiple backends.

    :param filename: The file you wish to open.
    :param full_header: Whether to read full headers.
    :param read_into_memory: Whether to read into memory.
    :param show_tqdm: Whether to display a progress bar.

    .. versionadded:: 1.0.2
    """
    if read_into_memory is None:
        read_into_memory = wc_c(filename) > 10 * 1024 * 1024
    if read_into_memory:
        return _MemoryAccessFastaView(
            filename=filename,
            full_header=full_header,
            show_tqdm=show_tqdm
        )
    else:
        return _DiskAccessFastaView(
            filename=filename,
            full_header=full_header,
            show_tqdm=show_tqdm
        )


def split_fasta(  # TODO: Add to commandline params
        fav: FastaViewType,
        out_dir_path: Optional[str] = None,
        safe_seqname: str = "convert"
):
    """
    Split input FASTA file into one-line FASTAs with one file per contig.

    :param fav: Source FASTA view.
    :param out_dir_path: Output directory.

    .. versionadded:: 1.0.2
    """
    out_dir_path = fav.filename + ".d" if out_dir_path is None else out_dir_path
    os.makedirs(out_dir_path, exist_ok=True)

    for seqname in fav.chr_names:
        safe_seqname = (
            seqname.
            replace(" ", "_").
            replace("\t", "_").
            replace("\\", "_").
            replace(":", "_").
            replace("*", "_").
            replace("?", "_").
            replace("\"", "_").
            replace("<", "_").
            replace(">", "_").
            replace("|", "_")
        )
        if seqname != safe_seqname:
            if safe_seqname == "convert":
                _lh.warning("seqname '%s' is not safe -- Converted to '%s'", seqname, safe_seqname)
            elif safe_seqname == "error":
                raise ValueError  # TODO
            elif safe_seqname == "skip":
                _lh.warning("seqname '%s' is not safe -- skipped", seqname)
                continue
        transcript_output_fasta = os.path.join(out_dir_path, f"{safe_seqname}.fa")
        with get_writer(transcript_output_fasta) as single_transcript_writer:
            single_transcript_writer.write(f">{seqname}\n{fav.sequence(seqname)}\n")
