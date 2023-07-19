"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

from abc import abstractmethod

from labw_utils.bioutils.datastructure.gv.gene import Gene
from labw_utils.typing_importer import Sequence


class GeneContainerInterface:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    @property
    @abstractmethod
    def number_of_genes(self) -> int:
        raise NotImplementedError

    @property
    @abstractmethod
    def gene_values(self) -> Sequence[Gene]:
        raise NotImplementedError

    @property
    @abstractmethod
    def gene_ids(self) -> Sequence[str]:
        raise NotImplementedError

    @abstractmethod
    def get_gene(self, gene_id: str) -> Sequence[Gene]:
        """
        Get gene by Gene ID.

        :return: Read-only Sequence of gene definitions with this gene ID.
            Will return empty sequence if not found.
        """
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
