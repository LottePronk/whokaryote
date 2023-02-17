#!/usr/bin/env python3
# Main script that wi160223ll be used on the commandline.

from whokaryote_scripts import *
import argparse
from pathlib import Path
import os
import time

parser = argparse.ArgumentParser(description="Classify metagenomic contigs as eukaryotic or prokaryotic")
parser.add_argument("--contigs", help="The path to your contigs file. It should be one (multi)fasta (DNA).")
parser.add_argument("--outdir", help="Specify the path to your preferred output directory. No / at the end.")
parser.add_argument("--prodigal_file", help=argparse.SUPPRESS)
parser.add_argument("--gff", help="If you already have gene predictions, specify path to the .gff file")
parser.add_argument("--f", action='store_true', help="If you want new multifastas with only eukaryotes and only "
                                                     "prokaryotes.")
parser.add_argument("--test", action='store_true', help="If you want to test it on a known dataset.")
parser.add_argument("--train", help="For training an RF on your own dataset. Provide name of RF output file.")
parser.add_argument("--minsize", default=5000, help="Select a minimum contig size in bp, default = 5000. Accuracy on \
contigs below 5000 is lower.")
#  parser.add_argument("--log", action='store_true', help="If you want a log file.")
parser.add_argument("--model", default="T", help="Choose the stand-alone model or the tiara-integrated model: S or T.\
 Option 'T' only works with argument --contigs")
parser.add_argument("--threads", default="1", help="Number of threads for Tiara to use.")
#  @TODO: integrate log file option into code..

args = parser.parse_args()

if args.prodigal_file:
    print("You are using the --prodigal_file option, which still works but is hidden and replaced with the --gff option.")
    args.gff = args.prodigal_file

if not args.outdir:
    print("Please specify an output directory with the option --outdir.")
    exit()

if os.path.isdir(args.outdir):
    print("Output directory is " + args.outdir)
else:
    try:
        os.mkdir(args.outdir)
    except PermissionError:
        print("There is no permission to change the directory. Please create output directory yourself.")

if args.contigs:
    #  filtered_contigs = 'empty'

    fasta_name = "contigs" + str(args.minsize) + ".fasta"
    filtered_contigs = os.path.join(args.outdir, fasta_name)

    if os.path.exists(filtered_contigs):
        print("Size-filtered contigs already exist. Using those.")
    else:
        try:
            print("Removing contigs with length <", args.minsize, "bp...")
            size_filter(args.contigs, args.outdir, size=int(args.minsize))
            fasta_name = "contigs" + str(args.minsize) + ".fasta"
            filtered_contigs = os.path.join(args.outdir, fasta_name)

        except NameError:
            print("Please provide contigs fasta file (nucleotide).")

    if args.model == "T":
        if os.path.exists((args.outdir + "/" + "tiara_pred.txt")):
            print("Tiara prediction file already present.")
        elif args.threads == "1":
            print("Model with tiara predictions selected.\nRunning tiara with " + args.threads + " thread...")
            run_tiara(filtered_contigs, args.outdir, args.threads)
        else:
            print("Model with tiara predictions selected.\nRunning tiara with " + args.threads + " threads...")
            run_tiara(filtered_contigs, args.outdir, args.threads)

    if args.gff:
        gene_predictions = args.gff
    else:
        print("Running prodigal...")
        prodigal_start = time.time()
        run_prodigal(filtered_contigs, args.outdir)
        print("Prodigal gene prediction successful. Saving gene coordinate file...")
        gene_predictions = os.path.join(args.outdir, "contigs_genes.gff")
        print("Gene coordinate file saved.")
        prodigal_tot = time.time() - prodigal_start
        print(f"Gene prediction took {prodigal_tot} seconds.")

if not args.contigs:
    print("No contig fasta was provided.")
    if args.gff:
        gene_predictions = args.gff
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

    print("Writing eukaryotic and prokaryotic contigs to separate fasta files.")
    filtered_contigs = os.path.join(args.outdir, "contigs" + str(args.minsize) + ".fasta")
    split_fasta_taxonomy(fasta_file=filtered_contigs, outdir=args.outdir)

    print("Writing contigs to separate fastas (eukaryotes.fasta, prokaryotes.fasta) successful. \nThese files were made"
          " from the size-filtered contigs and not from the original file.")
