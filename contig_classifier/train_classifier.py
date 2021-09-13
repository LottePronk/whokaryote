#!/usr/bin/env python3

""" A script that calculates the features of the contigs that are used in the classifier """
#  @TODO: Write the RF metrics to a 'metrics.txt' file.
#  @TODO: Generate unique identifiers for each training run (perhaps incl. date timestamp?)

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import seaborn as sn
import matplotlib.pyplot as plt
from sklearn import metrics
import joblib
from pathlib import Path
import os

contig_file = "unknown"


def add_tiara(dataframe, outdir):
    tiara_list = []

    tiara = os.path.join(outdir, "tiara_pred.txt")

    with open(tiara, newline='') as tiara_file:
        for line in tiara_file:

            if line.startswith("sequence_id"):
                continue

            all_line = line.strip().split("\t")

            seq_id = all_line[0].split(" ")[0]

            prediction = all_line[1]

            if prediction == "eukarya":
                prediction = "eukaryote"
            elif prediction == "archaea" or prediction == "bacteria" or prediction == "prokarya":
                prediction = "prokaryote"
            else:
                # . prediction == "organelle" or prediction == "unknown":
                prediction = prediction

            # . print(seq_id, prediction)

            tiara_list.append([seq_id, prediction])

    dataframe['tiara_pred'] = "other"

    for header in tiara_list:

        seq_id = header[0]

        if seq_id not in dataframe['contig'].values:
            pass

        else:
            predicted_class = header[1]
            dataframe.loc[dataframe['contig'] == seq_id, 'tiara_pred'] = predicted_class

    return dataframe


def calc_train_features(contig_file, outdir):
    # eukaryote_names = ["Phaseolus_vulgaris", "Puccinia_graminis", "Ustilago_maydis",
    #                    "Dictyostelium_discoideum", "Saprolegnia_parasitica",
    #                    "Amoebophrya_sp.", "Gigaspora_margarita", "Pyricularia_oryzae",
    #                    "Fusarium_oxysporum", "Aspergillus_oryzae", "Beauveria_bassiana",
    #                    "Phytophthora_parasitica", "Saccharomyces_cerevisiae"]
    #
    # prokaryote_names = ["Gemmata_obscuriglobus", "Streptomyces_albidoflavus",
    #                     "Candidatus_Prometheoarchaeum", "Clostridium_cellulovorans",
    #                     "Methanobacterium_subterraneum", "TPA_asm:_Burkholderia",
    #                     "Escherichia_coli", "Pseudomonas_protegens",
    #                     "Flavobacterium_lindanitolerans", "Chitinophaga_rhizosphaerae",
    #                     "Acidobacterium_capsulatum", "Methanococcus_maripaludis",
    #                     "Nostoc_punctiforme", "Xanthomonas_sacchari",
    #                     "Chloroflexus_aggregans", "Bacillus_thuringiensis"]

    data_dict = {"contig": [],
                 "organism": [],
                 "contig_length": [],
                 "nr_genes": [],
                 "kingdom": [],
                 "ratio_same_orientation": [],
                 "ID_general_avg": [],
                 "ID_general_std": [],
                 "gene_density": [],
                 "gene_length": [],
                 "ID_Q1": [],
                 "ID_Q3": []}

    gene_list = "empty"

    with open(contig_file, newline='') as coords_file:
        for line in coords_file:

            if line.startswith("# Seq"):

                if gene_list != "empty":  # This is the gene list from the previous contig
                    genes_same = 0  # <- <- or -> -> (++, --)
                    genes_diff = 0  # <- -> (-+)
                    genes_general = []  # All orientations together
                    data_dict["nr_genes"].append(len(gene_list))

                    if len(gene_list) < 2:
                        data_dict["ratio_same_orientation"].append(np.nan)

                        data_dict["ID_general_avg"].append(np.nan)
                        data_dict["ID_general_std"].append(np.nan)
                        data_dict["ID_Q1"].append(np.nan)
                        data_dict["ID_Q3"].append(np.nan)

                    if len(gene_list) >= 2:

                        for i in range(0, (len(gene_list) - 1)):
                            gene = gene_list[i]
                            next_gene = gene_list[i + 1]
                            distance = int(next_gene[1]) - int(gene[2])

                            if distance < 0:
                                distance = 0

                            gene_orientation = gene[3]
                            next_gene_orientation = next_gene[3]
                            orientation = str(gene_orientation) + "_" + str(next_gene_orientation)

                            genes_general.append(distance)

                            if orientation == "+_+" or orientation == "-_-":
                                genes_same += 1

                            if orientation == "-_+" or orientation == "+_-":
                                genes_diff += 1

                        ratio_same = genes_same / (genes_same + genes_diff)
                        data_dict["ratio_same_orientation"].append(ratio_same)

                        if len(genes_general) > 1:
                            data_dict["ID_general_avg"].append(np.mean(genes_general))
                            data_dict["ID_general_std"].append(np.std(genes_general))
                            data_dict["ID_Q1"].append(np.quantile(genes_general, 0.25))
                            data_dict["ID_Q3"].append(np.quantile(genes_general, 0.75))

                        if len(genes_general) == 1:
                            data_dict["ID_general_avg"].append(genes_general[0])
                            data_dict["ID_general_std"].append(0)
                            data_dict["ID_Q1"].append(0)
                            data_dict["ID_Q3"].append(0)

                        if len(genes_general) < 1:
                            data_dict["ID_general_avg"].append(np.nan)
                            data_dict["ID_general_std"].append(np.nan)
                            data_dict["ID_Q1"].append(np.nan)
                            data_dict["ID_Q3"].append(np.nan)

                    length = []
                    for gene in gene_list:
                        gene_length = int(gene[2]) - int(gene[1])
                        length.append(gene_length)

                    density = sum(length) / seqlength
                    data_dict["gene_density"].append(density)

                    if len(length) != 0:
                        avglength = sum(length) / len(length)
                        data_dict["gene_length"].append(avglength)

                    if len(length) == 0:
                        avglength = sum(length)
                        data_dict["gene_length"].append(avglength)

                header = line
                seqname = header.split(";")[2].split("=")[1].split(" ")[0]
                organism = "_".join(line.split(' ')[4:6])
                seqlength = int(header.split(";")[1].split("=")[1])
                domain = header.split('=')[3].split(' ')[0].split('_')[-1]

                data_dict["contig"].append(seqname)

                if domain == 'eukaryote':
                    data_dict["kingdom"].append(0)
                if domain == 'prokaryote':
                    data_dict["kingdom"].append(1)

                # if organism[0:(len(organism) - 2)] in eukaryote_names:
                #     data_dict["kingdom"].append(0)
                # if organism[0:(len(organism) - 2)] in prokaryote_names:
                #     data_dict["kingdom"].append(1)

                data_dict["organism"].append(organism)
                data_dict["contig_length"].append(seqlength)

                gene_list = []

            if line.startswith(">"):
                info = line.split("\n")[0].split(">")[1].split("_")
                gene_list.append(info)

        if gene_list != "empty":  # This is the gene list from the previous contig
            genes_same = 0  # <- <- or -> -> (++, --)
            genes_diff = 0  # <- -> (-+)
            genes_general = []  # All orientations together
            data_dict["nr_genes"].append(len(gene_list))

            if len(gene_list) < 2:
                data_dict["ratio_same_orientation"].append(np.nan)

                data_dict["ID_general_avg"].append(np.nan)
                data_dict["ID_general_std"].append(np.nan)
                data_dict["ID_Q1"].append(np.nan)
                data_dict["ID_Q3"].append(np.nan)

            if len(gene_list) >= 2:

                for i in range(0, (len(gene_list) - 1)):
                    gene = gene_list[i]
                    next_gene = gene_list[i + 1]
                    distance = int(next_gene[1]) - int(gene[2])

                    if distance < 0:
                        distance = 0

                    gene_orientation = gene[3]
                    next_gene_orientation = next_gene[3]
                    orientation = str(gene_orientation) + "_" + str(next_gene_orientation)

                    genes_general.append(distance)

                    if orientation == "+_+" or orientation == "-_-":
                        genes_same += 1

                    if orientation == "-_+" or orientation == "+_-":
                        genes_diff += 1

                ratio_same = genes_same / (genes_same + genes_diff)
                data_dict["ratio_same_orientation"].append(ratio_same)

                if len(genes_general) > 1:
                    data_dict["ID_general_avg"].append(np.mean(genes_general))
                    data_dict["ID_general_std"].append(np.std(genes_general))
                    data_dict["ID_Q1"].append(np.quantile(genes_general, 0.25))
                    data_dict["ID_Q3"].append(np.quantile(genes_general, 0.75))

                if len(genes_general) == 1:
                    data_dict["ID_general_avg"].append(genes_general[0])
                    data_dict["ID_general_std"].append(0)
                    data_dict["ID_Q1"].append(0)
                    data_dict["ID_Q3"].append(0)

                if len(genes_general) < 1:
                    data_dict["ID_general_avg"].append(np.nan)
                    data_dict["ID_general_std"].append(np.nan)
                    data_dict["ID_Q1"].append(np.nan)
                    data_dict["ID_Q3"].append(np.nan)

            length = []
            for gene in gene_list:
                gene_length = int(gene[2]) - int(gene[1])
                length.append(gene_length)

            density = sum(length) / seqlength
            data_dict["gene_density"].append(density)

            if len(length) != 0:
                avglength = sum(length) / len(length)
                data_dict["gene_length"].append(avglength)

            if len(length) == 0:
                avglength = sum(length)
                data_dict["gene_length"].append(avglength)

    df = pd.DataFrame(data_dict, columns=list(data_dict.keys()))

    df_tiara = add_tiara(df, outdir)

    original_nona = df_tiara.copy(deep=True)

    original_nona = original_nona.dropna()
    print("original nona", original_nona.shape)
    print("Check if dataframe is consistent...")
    for key in data_dict:
        print(key, len(data_dict[key]))

    del df['contig']
    del df['organism']
    del df['contig_length']
    del df['nr_genes']

    print("Check dataframe shape:", df.shape, "\n")

    features = pd.get_dummies(df)

    print("Print features.iloc[:,5:].head(5)")
    features.iloc[:, 5:].head(5)

    print("Describe features:")
    features.describe()

    features = features.dropna()
    print("Shape of features:", features.shape)

    # Labels are the values we want to predict
    labels = np.array(features['kingdom'])

    # Remove the labels from the features
    # axis 1 refers to the columns
    features = features.drop('kingdom', axis=1)

    # Convert to numpy array
    features = np.array(features)

    # Make random forest
    # Using Skicit-learn to split data into training and testing sets
    # Split the data into training and testing sets
    print("Splitting data into training and testing sets")
    train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size=0.30,
                                                                                random_state=42,
                                                                                stratify=labels)

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
    fig.figure.savefig(os.path.join(outdir, 'RF_confusionmatrix_5100500_proba_g3_3_tiarapred.png'))

    print('Accuracy: ', metrics.accuracy_score(test_labels, predictions))
    plt.show()

    print(metrics.classification_report(test_labels, predictions))

    print("Feature importances:")
    feature_importances = pd.Series(rf.feature_importances_).sort_values(ascending=False)
    print(feature_importances)

    #  sn.barplot(x=round(featureImportances,4), y=featureImportances)
    #  plt.xlabel('Features Importance')
    #  lt.savefig('feature_importance_5100500_proba.png')

    # Predict the whole dataset
    print("Predicting whole dataset")
    # all_predicted = rf.predict_proba(features)
    class_predicted = rf.predict(features)

    #  print(all_predicted)
    original_nona['predicted'] = class_predicted
    # df_new = pd.concat([original_nona, pd.DataFrame(class_predicted)], axis=1, join="inner")
    original_nona.to_csv(os.path.join(outdir, 'results_RF_5100500_g3_3_tiara.csv'))

    # save the RF model for later use
    joblib.dump(rf, os.path.join(Path(__file__).parents[1], "data", "random_forest5100500_g3_3_tiarapred.joblib"))
