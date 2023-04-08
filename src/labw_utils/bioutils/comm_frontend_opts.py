"""
frontend.py -- Utilities for other ``main`` functions.
"""

__all__ = (
    "FrontendOptSpec",
    "FrontendOptSpecs"
)

import argparse
from typing import Dict, Any, Mapping, Tuple

from labw_utils.devutils.decorators import copy_doc


class FrontendOptSpec:
    _args: Tuple[Any, ...]
    _kwargs: Mapping[str, Any]

    @copy_doc(argparse.ArgumentParser.add_argument)
    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs

    def patch(
            self,
            parser: argparse.ArgumentParser,
            **update_kwargs
    ) -> argparse.ArgumentParser:
        kwargs = dict(self._kwargs)
        kwargs.update(update_kwargs)
        parser.add_argument(*self._args, **kwargs)
        return parser

    @property
    def name(self) -> str:
        return self._args[0]


class FrontendOptSpecs:
    _inner_dict: Dict[str, FrontendOptSpec] = {}

    @staticmethod
    def add(opt_spec: FrontendOptSpec):
        FrontendOptSpecs._inner_dict[opt_spec.name] = opt_spec

    @staticmethod
    def patch(parser: argparse.ArgumentParser, name: str, **update_kwargs):
        return FrontendOptSpecs._inner_dict[name].patch(parser=parser, **update_kwargs)


FrontendOptSpecs.add(FrontendOptSpec(
    '-f',
    '--fasta',
    required=True,
    help="Path to input reference genome sequence, in FASTA format",
    nargs='?',
    type=str,
    action='store'
))

FrontendOptSpecs.add(FrontendOptSpec(
    '-g', '--gtf',
    required=True,
    help="Path to input genomic annotation in GTF format. Can be compressed.",
    nargs='?',
    type=str,
    action='store'
))