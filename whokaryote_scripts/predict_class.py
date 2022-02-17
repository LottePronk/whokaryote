import pandas as pd
import numpy as np
import joblib
from pathlib import Path
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


def predict_class(feature_path, outdir, model):

    if model == "T":
        print("Using model with tiara predictions to predict contig class...")
        #  model_file = "whokaryote_model_tiara.joblib"
        model_file = "RF_T_RBS_310122_bal.joblib"
    else:
        print("Using standard model to predict contig class...")
        #  model_file = "whokaryote_model_standard.joblib"
        model_file = "RF_S_RBS_240121_bal.joblib"

    feature_df = pd.read_csv(feature_path)

    if model == "T":
        df_tiara = add_tiara(feature_df, outdir)
        original_nona = df_tiara.copy(deep=True)
    else:
        original_nona = feature_df.copy(deep=True)

    original_nona = original_nona.dropna()

    del feature_df['contig']
    del feature_df['contig_length']
    del feature_df['nr_genes']

    #  features = pd.get_dummies(feature_df)
    features = feature_df.dropna()
    #  features = np.array(features)
    print("Used features:\n", features.describe())
    print(features.shape)
    print("Used model: ", model_file)
    loaded_rf = joblib.load(os.path.join(str(Path(__file__).parents[1]), "whokaryote_scripts/data", model_file))

    predictions = loaded_rf.predict(features)

    original_nona['predicted'] = predictions

    euk_namelist = original_nona[original_nona["predicted"] == 0]['contig'].to_list()
    prok_namelist = original_nona[original_nona["predicted"] == 1]['contig'].to_list()

    original_nona['predicted'] = original_nona['predicted'].replace([0, 1], ["eukaryote", "prokaryote"])

    csv_name = ("featuretable_predictions_" + model + ".tsv")
    original_nona.to_csv(os.path.join(outdir, csv_name), index=False, sep="\t")

    simple_table = original_nona[['contig', 'predicted']].copy()
    simple_output = ("whokaryote_predictions_" + model + ".tsv")
    simple_table.to_csv(os.path.join(outdir, simple_output), index=False, sep="\t")

    file = open(os.path.join(outdir, "eukaryote_contig_headers.txt"), 'w')
    for items in euk_namelist:
        file.writelines([items + '\n'])

    file.close()

    file = open(os.path.join(outdir, "prokaryote_contig_headers.txt"), 'w')
    for items in prok_namelist:
        file.writelines([items + '\n'])

    file.close()
