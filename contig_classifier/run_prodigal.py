""" Run gene prediction with prodigal, metagenomics mode. """

import subprocess


def run_prodigal(contig_file, outdir):
    genes_output = (outdir + "/contigs_genes.genes")
    proteins_output = (outdir + "/contigs_proteins.faa")

    subprocess.call(
        ["prodigal", "-i", contig_file, "-o", genes_output, "-a", proteins_output, "-p", "meta", "-f", "sco", "-q"])
