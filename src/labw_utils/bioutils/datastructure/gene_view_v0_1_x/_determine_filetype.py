"""
Common file suffixes. To be deprecated.
"""
import itertools

common_suffixes = {"GTF": (".gtf", ".gff"), "GFF3": (".gff3",)}

archive_suffixes = {
    "GZ": (".gz", ".gzip"),
    "LZMA": (".xz", ".lzma"),
    "LZ4": (".lz4",),
    "LZIP": (".lz",),
    "LZOP": (".lzop",),
    "BZ2": (".bz2",),
    "COMPRESS": (".z",),
    "BROTLI": (".brotli", ".br"),
    "ZSTD": (".zst", ".zstd"),
    "ZIP": (".zip",),
    "RAR": (".rar", ".rar5"),
    "7Z": (".7z",),
}


def get_file_type_from_suffix(filename: str) -> str:
    filename = filename.lower()
    archive_suffix_all_removed = False
    while not archive_suffix_all_removed:
        archive_suffix_all_removed = True
        for archive_suffix in itertools.chain(*archive_suffixes.values()):
            if filename.endswith(archive_suffix):
                filename = filename.rstrip(archive_suffix)
                archive_suffix_all_removed = False

    for file_type, real_suffixes in common_suffixes.items():
        for real_suffix in real_suffixes:
            if filename.endswith(real_suffix):
                return file_type
    return "UNKNOWN"
