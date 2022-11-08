import sys

import pysam


def extract(sam_filename: str, out_filename: str):
    written = 0
    with open(out_filename, "w") as writer:
        writer.write("\t".join(("chr", "s", "e")) + "\n")
        af = pysam.AlignmentFile(sam_filename, "rb")
        for aln_seg in af.fetch():
            writer.write("\t".join((
                aln_seg.reference_name,
                str(aln_seg.reference_start),
                str(aln_seg.reference_end)
            )) + "\n")
            written += 1
            if written % 100000 == 0:
                print(f"{written} records parsed")
        af.close()


if __name__ == "__main__":
    extract(sys.argv[1], sys.argv[2])
