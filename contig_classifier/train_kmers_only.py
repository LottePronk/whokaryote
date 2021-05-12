#!/usr/bin/env python3

""" A script that calculates the features of the contigs that are used in the classifier """

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import seaborn as sn
import matplotlib.pyplot as plt
from sklearn import metrics
import joblib
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
from skbio import Sequence
import itertools
import os
from pathlib import Path

contig_fasta = "/Users/lottepronk/Documents/Projects/metatoolkit/Data/510500/output_classifier/contigs5000.fasta"
outdir = "/Users/lottepronk/Documents/Projects/metatoolkit/Data/510500/kmer_classifier"
new_contig_fasta = "/Users/lottepronk/Documents/Projects/metatoolkit/Data/510500/output_classifier/contigs5000_upper.fasta"
contig_file = "/Users/lottepronk/Documents/Projects/metatoolkit/Data/510500/output_classifier/contigs_genes.genes"

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

eukaryote_names = ["Phaseolus_vulgaris", "Puccinia_graminis", "Ustilago_maydis",
                   "Dictyostelium_discoideum", "Saprolegnia_parasitica",
                   "Amoebophrya_sp.", "Gigaspora_margarita", "Pyricularia_oryzae",
                   "Fusarium_oxysporum", "Aspergillus_oryzae", "Beauveria_bassiana",
                   "Phytophthora_parasitica", "Saccharomyces_cerevisiae"]

prokaryote_names = ["Gemmata_obscuriglobus", "Streptomyces_albidoflavus",
                    "Candidatus_Prometheoarchaeum", "Clostridium_cellulovorans",
                    "Methanobacterium_subterraneum", "TPA_asm:_Burkholderia",
                    "Escherichia_coli", "Pseudomonas_protegens",
                    "Flavobacterium_lindanitolerans", "Chitinophaga_rhizosphaerae",
                    "Acidobacterium_capsulatum", "Methanococcus_maripaludis",
                    "Nostoc_punctiforme", "Xanthomonas_sacchari",
                    "Chloroflexus_aggregans", "Bacillus_thuringiensis"]

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
original_nona.to_csv(os.path.join(outdir, 'results_RF_510500_kmers.csv'))

original_nona = original_nona.dropna()
print("original nona", original_nona.shape)
print("Check if dataframe is consistent...")
for key in kmer_dict:
    print(key, len(kmer_dict[key]))

del kmer_df['contig']
del kmer_df['organism']

print("Check dataframe shape:", kmer_df.shape, "\n")

features = pd.get_dummies(kmer_df)

print("Print features.iloc[:,5:].head(5)")
features.iloc[:, 5:].head(5)

print("Describe features:")
features.describe()

features = features.dropna()
print("Shape of features:", features.shape)

original_table = features  # Going to use this for later

# Labels are the values we want to predict
labels = np.array(features['kingdom'])

# Remove the labels from the features
# axis 1 refers to the columns
features = features.drop('kingdom', axis=1)

# Saving feature names for later use
feature_list = list(features.columns)

# Convert to numpy array
features = np.array(features)

# Make random forest
# Using Skicit-learn to split data into training and testing sets
# Split the data into training and testing sets
print("Splitting data into training and testing sets")
train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size=0.30,
                                                                            random_state=42)

print('Training Features Shape:', train_features.shape)
print('Training Labels Shape:', train_labels.shape)
print('Testing Features Shape:', test_features.shape)
print('Testing Labels Shape:', test_labels.shape)

# Instantiate model with 100 decision trees
rf = RandomForestClassifier(n_estimators=100, random_state=42)

scores = cross_val_score(rf, train_features, train_labels, cv=5)
print("Cross-validation scores:", scores)
print("%0.2f accuracy with a standard deviation of %0.2f" % (scores.mean(), scores.std()))

# Train the model on training data
print("Training model")
rf.fit(train_features, train_labels)

# Use the forest's predict method on the test data
predictions = rf.predict(test_features)
# Calculate the absolute errors
errors = abs(predictions - test_labels)
# Print out the mean absolute error (mae)
print('Mean Absolute Error:', round(np.mean(errors), 2), 'degrees.')

confusion_matrix = pd.crosstab(test_labels, predictions, rownames=['Actual'], colnames=['Predicted'])
fig = sn.heatmap(confusion_matrix, annot=True)
fig.figure.savefig(os.path.join(outdir, 'RF_confusionmatrix_kmers_510500.png'))

print('Accuracy: ', metrics.accuracy_score(test_labels, predictions))
plt.show()

print(metrics.classification_report(test_labels, predictions))

print("Feature importances:")
featureImportances = pd.Series(rf.feature_importances_).sort_values(ascending=False)
print(featureImportances)

# Predict the whole dataset
print("Predicting whole dataset")
class_predicted = rf.predict(features)

#  print(all_predicted)
original_nona['predicted'] = class_predicted
# df_new = pd.concat([original_nona, pd.DataFrame(class_predicted)], axis=1, join="inner")
original_nona.to_csv(os.path.join(outdir, 'results_RF_kmers_510500.csv'))

# save the RF model for later use
joblib.dump(rf, os.path.join(Path(__file__).parents[1], "data", "random_forest_kmers_510500.joblib"))
