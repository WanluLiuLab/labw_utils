from __future__ import annotations

from abc import abstractmethod
from typing import Iterable

from labw_utils.bioutils.datastructure.gv.gene import Gene


class GeneContainerInterface:

    @property
    @abstractmethod
    def number_of_genes(self) -> int:
        raise NotImplementedError

    @property
    @abstractmethod
    def gene_values(self) -> Iterable[Gene]:
        raise NotImplementedError

    @property
    @abstractmethod
    def gene_ids(self) -> Iterable[str]:
        raise NotImplementedError

    @abstractmethod
    def get_gene(self, gene_id: str) -> Gene:
        raise NotImplementedError

    @abstractmethod
    def add_gene(self, gene: Gene) -> GeneContainerInterface:
        raise NotImplementedError

    @abstractmethod
    def del_gene(self, gene_id: str) -> GeneContainerInterface:
        raise NotImplementedError

    @abstractmethod
    def replace_gene(self, new_gene: Gene) -> GeneContainerInterface:
        raise NotImplementedError
