import argparse

from labw_utils.commonutils.serializer.toml import AbstractTOMLSerializable




class Reference(AbstractTOMLSerializable):
    ...


class RunSettings:
    _reference: Reference
    _bedtools_path: str
    _cache: str

    def run(self):
        self.create_cache_dir()

    def create_cache_dir(self):
        ...

    def merge_mask_bedfile(self):
        ...


    def mask_fasta(self):
        ...


if __name__ == '__main__':
    ...
