# Whokaryote

---
Classification of metagenomic contigs as eukaryotic/prokaryotic using biology-based features. 

Before predicting genes on your metagenomic contigs, use Whokaryote to determine which contigs need eukaryotic gene 
prediction and which need prokaryotic gene prediction.

Whokaryote uses a random forest classifier that uses gene-structure based features and Tiara
(https://github.com/ibe-uw/tiara) predictions to get highly accurate classification for most organisms. 
It shows improved accuracy on difficult genomes when compared to other classifiers.

---

## Installation

Whokaryote was developed to run from the commandline on a UNIX-based system such as Ubuntu or MacOS. 
If you want to use Whokaryote on your Windows PC, install Windows Subsystem for Linux (WSL). See
https://docs.microsoft.com/en-us/windows/wsl/about

###### Recommended installation:

We recommend that you use a fresh conda environment where you can install all other dependencies and tools.

1. Make a new and empty conda environment, and activate it:
- `conda create -n whokaryote`
- `conda activate whokaryote`
2. Install the following dependencies:
- `conda install -c bioconda prodigal`
- `pip install tiara` (see: https://github.com/ibe-uw/tiara)
3. Install whokaryote:
- Navigate to a directory you want to install whokaryote in.
- Clone whokaryote to this directory: `git clone https://git.wur.nl/lotte.pronk/whokaryote.git `
- `cd whokaryote`
- `python setup.py install`

---
## Using whokaryote

Use `whokaryote.py --help` to see the arguments you can use:
```
-h, --help            show this help message and exit
--contigs CONTIGS     The path to your contigs file. It should be one multifasta (DNA).
--outdir OUTDIR       Specify the path to your preferred output directory. No / at the end.
--prodigal_file PRODIGAL_FILE
                        If you already have prodigal gene predictions, specify path to the .genes or .gff file
--f                   If you want new multifastas with only eukaryotes and only prokaryotes. This can take a long time.
--test                If you want to test it on a known dataset.
--train               For training an RF on your own dataset
--minsize MINSIZE     Select a minimum contig size in bp, default = 5000. Accuracy oncontigs below 5000 is lower.
--model MODEL         Choose the stand-alone model or the tiara-integrated model: S or T. Option 'T' only works with argument --contigs
```

Warning! A tutorial on how to use --test and --train is still missing and will be added later.

### Example:

The most simple way to run Whokaryote with the tiara-integrated model.
If you already have gene prediction files (e.g. .gff from prodigal or .sco from prodigal), 
use the --prodigal_file option. If you don't have an annotation file yet, you can use only
--contigs and --outdir. Prodigal will then be run on your contigs file, which may take a while.

```
whokaryote.py --contigs contigs.fasta --outdir whokaryote_output --prodigal_file contigs_genes.gff
```

You will get the following output files:

- A file with all contig headers that were classified as eukaryotic (eukaryotic_contigs.txt), 
and a similar file for prokaryotic contig headers (prokaryotic_contigs.txt)
- A fasta file with only contigs that were longer than 5000 bp, called contigs5000.fasta
- A .tsv file with all the calculated features called calculated_features.tsv
- A file with the tiara predictions called tiara_pred.txt
