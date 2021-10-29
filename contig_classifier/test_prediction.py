#!/usr/bin/env python3

""" A script that is used to test the classifier. It can be used with a dataset for which it is known which contigs
 are eukaryotic and which ones are prokaryotic. """

import pandas as pd
import numpy as np
from pathlib import Path
import joblib
import os


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
                prediction = 0
            elif prediction == "archaea" or prediction == "bacteria" or prediction == "prokarya":
                prediction = 1
            else:
                # . prediction == "organelle" or prediction == "unknown":
                prediction = 2

            # . print(seq_id, prediction)

            tiara_list.append([seq_id, prediction])

    dataframe['tiara_pred'] = np.nan

    for header in tiara_list:

        seq_id = header[0]

        if seq_id not in dataframe['contig'].values:
            pass

        else:
            predicted_class = header[1]
            dataframe.loc[dataframe['contig'] == seq_id, 'tiara_pred'] = predicted_class

    return dataframe


def calc_test_features(contig_file, outdir, model):

    if model == "T":
        print("Using tiara-integrated model to predict contig class...")
        # model_file = "random_forest5100500_g3_3_tiarapred.joblib"
        model_file = "RF_T_291021.joblib"
    else:
        print("Using standard model to predict contig class...")
        model_file = "random_forest510500_g3_3.joblib"

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
                #  print(line)

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
                seqname = header.split(";")[2].split('"')[1].split(" ")[0]
                organism = "_".join(line.split(' ')[4:6])
                seqlength = int(header.split(";")[1].split("=")[1])
                domain = header.split('=')[3].split(' ')[0].split('_')[-1]

                data_dict["contig"].append(seqname)

                if domain == 'eukaryote':
                    data_dict["kingdom"].append(0)
                if domain == 'prokaryote':
                    data_dict["kingdom"].append(1)

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

    print("contig", len(data_dict["contig"]))
    print("organism", len(data_dict["organism"]))
    print("contig_length", len(data_dict["contig_length"]))
    print("nr_genes", len(data_dict["nr_genes"]))
    print("kingdom", len(data_dict["nr_genes"]))
    print("ratio_same_orientation", len(data_dict["ratio_same_orientation"]))
    print("ID_general_avg", len(data_dict["ID_general_avg"]))
    print("ID_general_std", len(data_dict["ID_general_std"]))
    print("gene_density", len(data_dict["gene_density"]))
    print("gene_length", len(data_dict["gene_length"]))
    print("ID_Q1", len(data_dict["ID_Q1"]))
    print("ID_Q3", len(data_dict["ID_Q3"]))

    df = pd.DataFrame(data_dict, columns=list(data_dict.keys()))

    df.to_csv(os.path.join(outdir, "df_features_noprediction.csv"), index=False)

    if model == "T":
        df_tiara = add_tiara(df, outdir)
        original_nona = df_tiara.copy(deep=True)
    else:
        original_nona = df.copy(deep=True)

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
    features = features.dropna()
    features = features.drop('kingdom', axis=1)
    features = np.array(features)
    print("Shape of features:", features.shape)

    loaded_rf = joblib.load(os.path.join(Path(__file__).parents[1], "data", model_file))

    predictions = loaded_rf.predict(features)

    original_nona['predicted'] = predictions

    import matplotlib.pyplot as plt
    from sklearn import metrics
    print('Accuracy: ', metrics.accuracy_score(original_nona["kingdom"], original_nona["predicted"]))
    plt.show()
    #  @TODO: write this information to a metrics.txt file

    print(metrics.classification_report(original_nona["kingdom"], original_nona["predicted"]))

    euk_namelist = original_nona[original_nona["predicted"] == 0]['contig'].to_list()
    prok_namelist = original_nona[original_nona["predicted"] == 1]['contig'].to_list()

    original_nona.to_csv(os.path.join(outdir, "test_predictions.csv"), index=False)

    file = open(os.path.join(outdir, "eukaryote_contig_headers.txt"), 'w')
    for items in euk_namelist:
        file.writelines([items + '\n'])

    file.close()

    file = open(os.path.join(outdir, "prokaryote_contig_headers.txt"), 'w')
    for items in prok_namelist:
        file.writelines([items + '\n'])

    file.close()
