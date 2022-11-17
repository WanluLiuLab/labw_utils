from labw_utils.bioutils.datastructure.gene_tree import GeneTree
from labw_utils.bioutils.parser.gtf import GtfIterator, GtfIteratorWriter


def test():
    GtfIteratorWriter.write_iterator(
        GeneTree.from_feature_iterator(
            GtfIterator("/home/yuzj/Desktop/BioRef/chr1.ncbiRefSeq_h100000.gtf")
        ).to_feature_iterator(),
        filename="/home/yuzj/Desktop/BioRef/chr1.ncbiRefSeq_h100000.std.gtf"
    )


if __name__ == "__main__":
    test()
