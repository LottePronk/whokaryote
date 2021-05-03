#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 2021

@author: lottepronk
"""

#  from calculate_features import calc_features
import pandas as pd
import numpy as np
import joblib
from pathlib import Path


def predict_class(feature_path, outdir):

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

    loaded_rf = joblib.load(str(Path(__file__).parents[1] / "data/random_forest.joblib"))

    predictions = loaded_rf.predict(features)

    original_nona['predicted'] = predictions

    euk_namelist = original_nona[original_nona["predicted"] == 0]['contig'].to_list()
    prok_namelist = original_nona[original_nona["predicted"] == 1]['contig'].to_list()

    #  original_nona.to_csv(str(Path(__file__).parents[1] / "data/featuretable_predictions.csv"), index=False)

    original_nona.to_csv(outdir + "/featuretable_predictions.csv", index=False)

    file = open((outdir + "/eukaryote_contig_headers.txt"), 'w')
    for items in euk_namelist:
        file.writelines([items + '\n'])

    file.close()

    file = open((outdir + "/prokaryote_contig_headers.txt"), 'w')
    for items in prok_namelist:
        file.writelines([items + '\n'])

    file.close()
