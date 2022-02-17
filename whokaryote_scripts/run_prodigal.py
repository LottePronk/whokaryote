""" Run gene prediction with prodigal, metagenomics mode. """

import subprocess
import os


def run_prodigal(contig_file, outdir):
    genes_output = os.path.join(outdir, "contigs_genes.genes")
    proteins_output = os.path.join(outdir, "contigs_proteins.faa")

    subprocess.call(
        ["prodigal", "-i", contig_file, "-o", genes_output, "-a", proteins_output, "-p", "meta", "-f", "gff", "-q"])
