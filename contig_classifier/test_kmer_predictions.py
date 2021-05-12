import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
import joblib
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
from skbio import Sequence
import itertools
import os


eukaryote_names = ["Arabidopsis_thaliana", "Candida_dubliniensis", "Drosophila_melanogaster", "Homo_sapiens",
                   "Neurospora_crassa", "Nosema_ceranae", "Physcomitrium_patens", "Synchytrium_microbalum",
                   "Toxoplasma_gondii", "Tremella_mesenterica", "TPA_asm:_Arabidopsis", "Physcomitrella_patens"]

prokaryote_names = ["Bacteroides_heparinolyticus", "Chlamydia_pecorum", "Deinococcus_radiodurans",
                    "Gemmatimonas_phototrophica", "Myxococcus_xanthus", "Nocardia_mangyaensis",
                    "Rickettsia_heilongjiangensis", "Sphingobacterium_lactis"]

contig_fasta = "/Users/lottepronk/Documents/Projects/metatoolkit/Data/new_test/g2contigs5100500.fasta"
outdir = "/Users/lottepronk/Documents/Projects/metatoolkit/Data/new_test/output_kmer"
new_contig_fasta = "/Users/lottepronk/Documents/Projects/metatoolkit/Data/new_test/g2contigs5100500_upper.fasta"
contig_file = "/Users/lottepronk/Documents/Projects/metatoolkit/Data/new_test/contigs_genes.genes"

records = (rec.upper() for rec in SeqIO.parse(contig_fasta, "fasta"))
count = SeqIO.write(records, new_contig_fasta, "fasta")
print("Converted %i records to upper case" % count)

kmer_dict = {"contig": [],
             "organism": [],
             "kingdom": []}

bases = ['A', 'T', 'G', 'C', 'N']
k = 3
kmer_list = [''.join(p) for p in itertools.product(bases, repeat=k)]

for triplet in kmer_list:
    kmer_dict[triplet] = []

not_base = ['U', '(i)', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', '-']

with open(new_contig_fasta, newline='') as contigs:

    for title, seq in SimpleFastaParser(contigs):
        contig = title.split(' ')[0]
        organism = title.split(' ')[1]
        kmer_dict["contig"].append(contig)
        kmer_dict["organism"].append(organism)

        if organism in eukaryote_names:
            kmer_dict["kingdom"].append(0)
        if organism in prokaryote_names:
            kmer_dict["kingdom"].append(1)

        for char in not_base:
            seq = seq.replace(char, 'N')
        s = Sequence(seq)
        freqs = s.kmer_frequencies(3, relative=True, overlap=False)

        for kmer in kmer_list:
            if kmer not in freqs:
                freqs[kmer] = 0.00

        for triplet, freq in freqs.items():
            kmer_dict[str(triplet)].append(freq)

kmer_df = pd.DataFrame(kmer_dict, columns=list(kmer_dict.keys()))

###########

###########
original_nona = kmer_df.copy(deep=True)
original_nona.to_csv(os.path.join(outdir, 'results_g2_510500_kmers.csv'))

original_nona = original_nona.dropna()
print("original nona", original_nona.shape)
print("Check if dataframe is consistent...")
for key in kmer_dict:
    print(key, len(kmer_dict[key]))

del kmer_df['contig']
del kmer_df['organism']

print("Check dataframe shape:", kmer_df.shape, "\n")

features = pd.get_dummies(kmer_df)
features = features.dropna()
features = features.drop('kingdom', axis=1)
features = np.array(features)
print("Shape of features:", features.shape)

loaded_rf = joblib.load(os.path.join(Path(__file__).parents[1], "data", "random_forest_kmers_510500.joblib"))

predictions = loaded_rf.predict(features)

original_nona['predicted'] = predictions

print('Accuracy: ', metrics.accuracy_score(original_nona["kingdom"], original_nona["predicted"]))
plt.show()

print(metrics.classification_report(original_nona["kingdom"], original_nona["predicted"]))

euk_namelist = original_nona[original_nona["predicted"] == 0]['contig'].to_list()
prok_namelist = original_nona[original_nona["predicted"] == 1]['contig'].to_list()

original_nona.to_csv(os.path.join(outdir, "test_predictions_kmers.csv"), index=False)
