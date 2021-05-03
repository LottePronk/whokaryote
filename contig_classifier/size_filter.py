import Bio
from Bio.SeqIO.FastaIO import SimpleFastaParser


def size_filter(contig_file, outdir, size):
    with open(contig_file) as in_handle:
        with open((outdir + "/contigs5000.fasta"), 'w', newline='') as filtered_contigs:
            for title, seq in SimpleFastaParser(in_handle):
                if len(seq) >= size:
                    filtered_contigs.write(">%s\n%s\n" % (title, seq))