from contig_classifier import *
import argparse
from pathlib import Path
import os
import time

parser = argparse.ArgumentParser("Classify metagenomic contigs as eukaryotic or prokaryotic")
parser.add_argument("--contigs", help="The path to your contigs file. It should be one multifasta.")
parser.add_argument("--outdir", help="Specify the path to your preferred output directory. No / at the end.")
parser.add_argument("--prodigal_file", help="If you already have prodigal gene predictions, specify path to the "
                                            ".genes or .gff file")
parser.add_argument("--f", action='store_true', help="If you want new multifastas with only eukaryotes and only "
                                                     "prokaryotes. This can take a long time.")
parser.add_argument("--test", action='store_true', help="If you want to test it on a known dataset.")
parser.add_argument("--train", action='store_true', help="For training an RF on your own dataset")
parser.add_argument("--minsize", default=5000, help="Select a minimum contig size in bp, default = 5000. Accuracy on\
contigs below 5000 is lower.")
parser.add_argument("--log", action='store_true', help="If you want a log file.")
#  @TODO: integrate log file option into code

args = parser.parse_args()

if args.contigs:
    print("Removing contigs with length <", args.minsize, "bp...")
    size_filter(args.contigs, args.outdir, size=args.minsize)
    fasta_name = "contigs" + str(args.minsize) + ".fasta"

    filtered_contigs = os.path.join(args.outdir, fasta_name)

    run_tiara(filtered_contigs, args.outdir)

    if not args.prodigal_file:
        print("Running prodigal...")
        run_prodigal(filtered_contigs, args.outdir)
        print("Prodigal successful. Saving gene coordinate file...")
        contig_file = os.path.join(args.outdir, "contigs_genes.genes")
        print("Gene coordinate file saved.")

        if not args.test and not args.train:
            print("Calculating features...")
            calc_features(contig_file, args.outdir)
            print("Calculating features successful.")

        if args.test:
            print("Calculating features and predicting test data...")
            calc_test_features(contig_file, args.outdir)
            print("Test successful...")

        if args.train:
            print("Training a new classifier...")
            calc_train_features(contig_file, args.outdir)
            print("Training successful...")

if args.prodigal_file:

    print("Calculating features from gene coordinates file...")
    if not args.test and not args.train:
        calc_features(args.prodigal_file, args.outdir)
        print("Calculating features successful.")

    if args.test:
        print("Testing classifier on dataset with known taxonomy...")
        calc_test_features(args.prodigal_file, args.outdir)
        print("Testing successful. Check the output files.")

    if args.train:
        print("Training a new classifier...")
        calc_train_features(args.prodigal_file, args.outdir)
        print("Training successful...")

if not args.test and not args.train:
    print("Predicting contig class...")
    feature_path = os.path.join(args.outdir, "featuretable.csv")
    predict_class(feature_path, args.outdir)
    print("Prediction successful! See output directory.")

if args.f:
    print("Writing eukaryotic and prokaryotic contigs to separate fasta files. This can take very long...")
    script_path = os.path.join(str(Path(__file__).parents[1]), "contig_classifier", "get_euk_prok_fasta.sh")
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
