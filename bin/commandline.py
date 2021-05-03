from contig_classifier import *
import argparse
from pathlib import Path

parser = argparse.ArgumentParser("Classify metagenomic contigs as eukaryotic or prokaryotic")
parser.add_argument("--contigs", help="The path to your contigs file. It should be one multifasta.")
parser.add_argument("--outdir", help="Specify the path to your preferred output directory. No / at the end.")
parser.add_argument("--prodigal_file", help="If you already have prodigal gene predictions, specify path to the "
                                            ".genes or .gff file")
parser.add_argument("--f", action='store_true', help="If you want new multifastas with only eukaryotes and only "
                                                     "prokaryotes. This can take a long time.")
parser.add_argument("--test", action='store_true', help="If you want to test it on a known dataset.")
parser.add_argument("--train", action='store_true', help="For training an RF on your own dataset")

args = parser.parse_args()

print("Removing contigs with length < 5000 bp...")
size_filter(args.contigs, args.outdir, size=5000)

filtered_contigs = str(args.outdir + "/contigs5000.fasta")

if not args.prodigal_file:
    print("Running prodigal...")
    run_prodigal(filtered_contigs, args.outdir)
    print("Prodigal successful. Saving gene coordinate file...")
    contig_file = (args.outdir + "/contigs_genes.genes")
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
    feature_path = (args.outdir + "/featuretable.csv")
    predict_class(feature_path, args.outdir)
    print("Prediction successful! See output directory.")

if args.f:
    print("Writing eukaryotic and prokaryotic contigs to separate fasta files. This can take very long...")
    script_path = str(Path(__file__).parents[1] / "contig_classifier/get_euk_prok_fasta.sh")
    input_file = args.contigs
    output_file = (args.outdir + "/lin_contigs.fasta")
    euk_headers = (args.outdir + "/eukaryote_contig_headers.txt")
    prok_headers = (args.outdir + "/prokaryote_contig_headers.txt")
    euk_fasta = args.outdir + "/eukaryotic_contigs.fasta"
    prok_fasta = args.outdir + "/prokaryotic_contigs.fasta"

    subprocess.run(
        ["bash", script_path,
         input_file,
         output_file,
         euk_headers,
         prok_headers,
         euk_fasta,
         prok_fasta])

    print("Writing contigs to separate fastas (eukaryote, prokaryote) successful.")
