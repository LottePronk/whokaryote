# Whokaryote

Classification of metagenomic contigs as eukaryotic/prokaryotic using biology-based features.

Whokaryote uses a random forest classifier that uses gene-structure based features and Tiara
(https://github.com/ibe-uw/tiara) predictions to predict whether a contig is from a eukaryote or from a prokaryote.

You can use Whokaryote to determine which contigs need eukaryotic gene prediction and which need prokaryotic gene prediction.

---

## Installation

Whokaryote was developed to run from the commandline on a UNIX-based system such as Ubuntu or MacOS. 
If you want to use Whokaryote on your Windows PC, install Windows Subsystem for Linux (WSL). See
https://docs.microsoft.com/en-us/windows/wsl/about

###### Recommended installation:

We recommend that you install whokaryote and its dependencies in a new conda environment.
You can download miniconda here: https://docs.conda.io/en/latest/miniconda.html

1. Make a new and empty conda environment, and activate it:
- `conda create -n whokaryote python==3.8`
- `conda activate whokaryote`
2. Install dependencies with these commands:
- `conda install -c bioconda prodigal`
- `conda install pip`
- `python -m pip install tiara` (see: https://github.com/ibe-uw/tiara)
3. Install whokaryote:
- Navigate to a directory you want to install whokaryote in.
- Clone whokaryote to this directory: `git clone https://git.wur.nl/lotte.pronk/whokaryote.git `
- Go to the whokaryote directory: `cd whokaryote`
- Install whokaryote in your (conda) environment with`python setup.py install`

---
## Using whokaryote

Use `whokaryote.py --help` to see all the options:
```
-h, --help            show this help message and exit
--contigs CONTIGS     The path to your contigs file. It should be one multifasta (DNA).
--outdir OUTDIR       Specify the path to your preferred output directory. No / at the end.
--prodigal_file PRODIGAL_FILE
                        If you already have prodigal gene predictions, specify path to the .genes or .gff file
--f                   If you want new multifastas with only eukaryotes and only prokaryotes. This can take a long time.
--minsize MINSIZE     Select a minimum contig size in bp, default = 5000. Accuracy oncontigs below 5000 is lower.
--model MODEL         Choose the stand-alone model or the tiara-integrated model: S or T. Option 'T' only works with argument --contigs
```

### Example:

```
whokaryote.py --contigs contigs.fasta --outdir whokaryote_output --prodigal_file contigs_genes.gff
```
This is the standard way to run Whokaryote. If you don't specify --model, it will run the model that integrates 
Tiara predictions. This is recommended. 

If you already have gene annotation files (.gff or gene coordinates file from prodigal), 
use the --prodigal_file option. 

If you don't have an annotation file yet, you can use only
--contigs and --outdir. Prodigal will then be run on your contigs file, which may take a while.

With the example command, you will get the following output files:

- A file with all contig headers that were classified as eukaryotic (eukaryotic_contigs.txt), 
and a similar file for prokaryotic contig headers (prokaryotic_contigs.txt)
- A fasta file with only contigs that were longer than 5000 bp, called contigs5000.fasta
- A .tsv file with all the calculated features called calculated_features.tsv
- A file with the tiara predictions called tiara_pred.txt
