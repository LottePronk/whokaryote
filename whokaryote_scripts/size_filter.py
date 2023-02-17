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
                    print("Warning! sequence" + title + "is longer than 32000000 BP. If you did not provide a GFF"
                     "file, Prodigal will be used to predict genes. Prodigal does not accept sequences > 32 mbp!")

            print("Total contigs checked:\t" + str(total_contigs),
                  '\nNumber of contigs >=' + str(size) + 'bp:\t' + str(kept_contigs),
                  '\nNumber of contigs <' + str(size) + 'bp and removed:\t' + str(removed_contigs))

            if kept_contigs == 0:
                sys.exit("No contigs longer than" + str(size) + "bp detected.\nYour contigs cannot be classified.")

def split_fasta_taxonomy(fasta_file, outdir):
    """
    Split fasta file into separate files for each taxonomy.
    output: outdir/eukaryote.fasta, outdir/prokaryote.fasta
    """
    eukaryote_headers = os.path.join(outdir, "eukaryote_contig_headers.txt")
    prokaryote_headers = os.path.join(outdir, "prokaryote_contig_headers.txt")

    euk_file = open(eukaryote_headers, "r")
    euk_headers = euk_file.read()
    euk_headers_list = euk_headers.split("\n")  # list of all eukaryote headers
    euk_file.close()

    prok_file = open(prokaryote_headers, "r")
    prok_headers = prok_file.read()
    prok_headers_list = prok_headers.split("\n")  # list of all prokaryote headers
    prok_file.close()

    with open(fasta_file, "r") as f:
        for title, seq in SimpleFastaParser(f):
            title_id = title.split(" ")[0]
            if title_id in euk_headers_list:
                with open(os.path.join(outdir, "eukaryotes.fasta"), "a") as euk:
                    euk.write(">" + title + "\n" + seq + "\n")
            elif title_id in prok_headers_list:
                with open(os.path.join(outdir, "prokaryotes.fasta"), "a") as prok:
                    prok.write(">" + title + "\n" + seq + "\n")
            else:
                with open(os.path.join(outdir, "unclassified.fasta"), "a") as unclass:
                    unclass.write(">" + title + "\n" + seq + "\n")
