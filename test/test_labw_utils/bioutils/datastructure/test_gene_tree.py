from labw_utils.bioutils.datastructure.gene_tree import GeneTree
from labw_utils.bioutils.parser.gtf import GtfIterator


def test():
    gt = GeneTree.from_feature_iterator(
        GtfIterator("/home/yuzj/Desktop/BioRef/chr1.ncbiRefSeq.gtf")
    )


if __name__ == "__main__":
    test()
