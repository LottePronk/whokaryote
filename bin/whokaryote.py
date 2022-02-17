# Main script that will be used on the commandline.

from whokaryote_scripts import *
import argparse
from pathlib import Path
import os
import time

parser = argparse.ArgumentParser("Classify metagenomic contigs as eukaryotic or prokaryotic")
parser.add_argument("--contigs", help="The path to your contigs file. It should be one multifasta (DNA).")
parser.add_argument("--outdir", help="Specify the path to your preferred output directory. No / at the end.")
parser.add_argument("--prodigal_file", help="If you already have prodigal gene predictions, specify path to the "
                                            ".genes or .gff file")
parser.add_argument("--f", action='store_true', help="If you want new multifastas with only eukaryotes and only "
                                                     "prokaryotes. This can take a long time.")
parser.add_argument("--test", action='store_true', help="If you want to test it on a known dataset.")
parser.add_argument("--train", help="For training an RF on your own dataset. Provide name of RF output file.")
parser.add_argument("--minsize", default=5000, help="Select a minimum contig size in bp, default = 5000. Accuracy on\
contigs below 5000 is lower.")
#  parser.add_argument("--log", action='store_true', help="If you want a log file.")
parser.add_argument("--model", default="T", help="Choose the stand-alone model or the tiara-integrated model: S or T.\
 Option 'T' only works with argument --contigs")
parser.add_argument("--threads", default="1", help="Number of threads for Tiara to use.")
#  @TODO: integrate log file option into code.
#  @TODO: Check why multithreading with Tiara makes the classification slower.

args = parser.parse_args()

if os.path.isdir(args.outdir):
    print("Output directory is " + args.outdir)
else:
    try:
        os.mkdir(args.outdir)
    except PermissionError:
        print("There is no permission to change the directory. Please create output directory yourself.")

if args.contigs:
    filtered_contigs = 'empty'
    try:
        print("Removing contigs with length <", args.minsize, "bp...")
        size_filter(args.contigs, args.outdir, size=int(args.minsize))
        fasta_name = "contigs" + str(args.minsize) + ".fasta"

        filtered_contigs = os.path.join(args.outdir, fasta_name)
    except NameError:
        print("Please provide contigs fasta file (nucleotide).")

    if args.model == "T":
        if args.threads == "1":
            print("Model with tiara predictions selected.\nRunning tiara with " + args.threads + " thread...")
        else:
            print("Model with tiara predictions selected.\nRunning tiara with " + args.threads + " threads...")
        run_tiara(filtered_contigs, args.outdir, args.threads)

    if args.prodigal_file:
        gene_predictions = args.prodigal_file
    else:
        print("Running prodigal...")
        prodigal_start = time.time()
        run_prodigal(filtered_contigs, args.outdir)
        print("Prodigal successful. Saving gene coordinate file...")
        gene_predictions = os.path.join(args.outdir, "contigs_genes.genes")
        print("Gene coordinate file saved.")
        prodigal_tot = time.time() - prodigal_start
        print(f"Gene prediction took {prodigal_tot} seconds.")

if not args.contigs:
    print("No contig fasta was provided.")
    if args.prodigal_file:
        gene_predictions = args.prodigal_file
    else:
        print("No gene predictions found. Provide contig fasta or prodigal output file.")
        quit()

if not args.test and not args.train:
    print("Calculating gene structure features with gene predictions...")
    calc_start = time.time()
    calc_features(gene_predictions, args.outdir)
    print("Calculating features successful.")
    calc_tot = time.time() - calc_start
    print(f"Calculating gene structure features took {calc_tot} seconds.")

    print("Predicting contig class...")
    predict_start = time.time()
    feature_path = os.path.join(args.outdir, "featuretable.csv")
    predict_class(feature_path, args.outdir, args.model)
    print("Prediction successful! See output directory.")
    predict_tot = time.time() - calc_start
    print(f"Predicting contig class took {predict_tot} seconds.")

if args.test:
    pass
    # print("Calculating features and predicting test data...")
    # test_start = time.time()
    # calc_test_features(gene_predictions, args.outdir, args.model)
    # print("Test successful...")
    # test_tot = time.time() - test_start
    # print(f"Prediction of test data took {test_tot} seconds.")

if args.train:
    pass
    print("Training a new classifier...")
    calc_training_features(gene_predictions, args.outdir, args.train)
    print("Training successful...")

if args.f:
    print("Writing eukaryotic and prokaryotic contigs to separate fasta files. This can take very long...")
    script_path = os.path.join(str(Path(__file__).parents[1]), "whokaryote_scripts", "get_euk_prok_fasta.sh")
    input_file = args.contigs
    output_file = os.path.join(args.outdir, "lin_contigs.fasta")
    euk_headers = os.path.join(args.outdir, "eukaryote_contig_headers.txt")
    prok_headers = os.path.join(args.outdir, "prokaryote_contig_headers.txt")
    euk_fasta = os.path.join(args.outdir, "eukaryotic_contigs.fasta")
    prok_fasta = os.path.join(args.outdir, "prokaryotic_contigs.fasta")

    subprocess.run(
        ["bash", script_path,
         input_file,
         output_file,
         euk_headers,
         prok_headers,
         euk_fasta,
         prok_fasta])

    print("Writing contigs to separate fastas (eukaryote, prokaryote) successful.")
