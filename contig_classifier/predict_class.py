import pandas as pd
import numpy as np
import joblib
from pathlib import Path
import os


def predict_class(feature_path, outdir, model):

    if model == "T":
        print("Using model with tiara predictions to predict contig class...")
        model_file = "random_forest5100500_g3_3_tiarapred.joblib"
    else:
        print("Using standard model to predict contig class...")
        model_file = "random_forest510500_g3_3.joblib"

    feature_df = pd.read_csv(feature_path)

    original_nona = feature_df.copy(deep=True)
    original_nona = original_nona.dropna()

    del feature_df['contig']
    del feature_df['contig_length']
    del feature_df['nr_genes']

    features = pd.get_dummies(feature_df)
    features = features.dropna()
    features = np.array(features)
    print(features.shape)

    loaded_rf = joblib.load(os.path.join(Path(__file__).parents[1], "data", model_file))

    predictions = loaded_rf.predict(features)

    original_nona['predicted'] = predictions

    euk_namelist = original_nona[original_nona["predicted"] == 0]['contig'].to_list()
    prok_namelist = original_nona[original_nona["predicted"] == 1]['contig'].to_list()

    csv_name = ("featuretable_predictions_" + model + ".csv")
    original_nona.to_csv(os.path.join(outdir, csv_name), index=False)

    file = open(os.path.join(outdir, "eukaryote_contig_headers.txt"), 'w')
    for items in euk_namelist:
        file.writelines([items + '\n'])

    file.close()

    file = open(os.path.join(outdir, "prokaryote_contig_headers.txt"), 'w')
    for items in prok_namelist:
        file.writelines([items + '\n'])

    file.close()
