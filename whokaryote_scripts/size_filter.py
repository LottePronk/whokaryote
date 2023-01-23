import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os


def size_filter(contig_file, outdir, size):
    file_name = "contigs" + str(size) + ".fasta"
    with open(contig_file) as in_handle:
        with open((os.path.join(outdir, file_name)), 'w', newline='') as filtered_contigs:
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
                if len(seq) > 32000000:
                    print("Warning! sequence" + title + "is longer than 32000000 BP. If you did not provide a GFF\
                     file, Prodigal will be used to predict genes. Prodigal does not accept sequences > 32 mbp!")

            print("Total contigs checked:\t" + str(total_contigs),
                  '\nNumber of contigs >=' + str(size) + 'bp:\t' + str(kept_contigs),
                  '\nNumber of contigs <' + str(size) + 'bp and removed:\t' + str(removed_contigs))

            if kept_contigs == 0:
                sys.exit("No contigs longer than" + str(size) + "bp detected.\nYour contigs cannot be classified.")

def split_fasta_taxonomy(fasta_file, outdir, headerfile):
    """
    Split fasta file into separate files for each taxonomy.
    output: outdir/eukaryote.fasta, outdir/prokaryote.fasta
    """

    euk_file = open(headerfile, "r")
    euk_headers = euk_file.read()
    euk_headers_list = euk_headers.split("\n")  # list of all eukaryote headers
    euk_file.close()

    with open(fasta_file, "r") as f:
        for title, seq in SimpleFastaParser(f):
            if title in euk_headers_list:
                with open(os.path.join(outdir, "eukaryotes.fasta"), "a") as g:
                    g.write(">" + title + "\n" + seq + "\n")
            else:
                with open(os.path.join(outdir, "prokaryotes.fasta"), "a") as g:
                    g.write(">" + title + "\n" + seq + "\n")