import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os


def size_filter(contig_file, outdir, size):
    with open(contig_file) as in_handle:
        with open((os.path.join(outdir, "contigs5000.fasta")), 'w', newline='') as filtered_contigs:
            total_contigs = 0
            kept_contigs = 0
            removed_contigs = 0

            for title, seq in SimpleFastaParser(in_handle):
                total_contigs += 1

                if len(seq) >= size:
                    filtered_contigs.write(">%s\n%s\n" % (title, seq))
                    kept_contigs += 1
                if len(seq) < size:
                    removed_contigs += 1

            print("Total contigs checked:\t" + str(total_contigs),
                  '\nNumber of contigs >= 5000 bp:\t' + str(kept_contigs),
                  '\nNumber of contigs < 5000 bp and removed:\t' + str(removed_contigs))

            if kept_contigs == 0:
                sys.exit("No contigs longer than 5000 bp detected.\nYour contigs cannot be classified.")
