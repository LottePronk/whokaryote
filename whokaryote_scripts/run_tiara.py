""" Run Tiara, default settings. """

import subprocess
import os


def run_tiara(contig_file, outdir, threads):
    prediction_output = os.path.join(outdir, "tiara_pred.txt")

    subprocess.call(
        ["tiara", "-i", contig_file, "-o", prediction_output, "-m 1000", "-t", threads])
