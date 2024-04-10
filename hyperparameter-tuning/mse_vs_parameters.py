#!/usr/bin/env python3
import pandas as pd
import numpy as np
import json
import os
import argparse

import tensorflow as tf
import tensorflow.keras as keras
from keras import layers


def TrainModel(pdf, hidden_nodes):
    print("Training Model")
    features = ["jet_pt_raw","jet_nparts","jet_area","jet_angularity","jet_track_pt_0",
                        "jet_track_pt_1","jet_track_pt_2","jet_track_pt_3","jet_track_pt_4",
                        "jet_track_pt_5","jet_track_pt_6","jet_track_pt_7"]
    target = "jet_pt_truth" 
    if len(hidden_nodes) != 3:
       print("Error: hidden_nodes must be a list of length 3")
       return
   
    inputlen = len(features)
    model = keras.Sequential(
            [
                layers.Dense(hidden_nodes[0],input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(hidden_nodes[1], activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(hidden_nodes[2], activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(1, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
            ])
    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    print(f"Fitting Model: {hidden_nodes}")
    model.fit(pdf[features], pdf[target], epochs=25, batch_size=128, verbose=0, validation_split=0.2)
    
    model_mse = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
    model_trainable_params = 0
    model_trainable_params += hidden_nodes[0]*(inputlen+1)
    model_trainable_params += hidden_nodes[1]*(hidden_nodes[0]+1)
    model_trainable_params += hidden_nodes[2]*(hidden_nodes[1]+1)
    model_trainable_params += 1*(hidden_nodes[2]+1)
    
    return model_mse, model_trainable_params
    
def GetPDF():
    print("Getting PDF")
    file = "/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/root-files/AuAu_R04/test/AuAu_R04_50_test_sample.h5"
    pdf = pd.read_hdf(file, key='df')
    features = ["jet_pt_raw","jet_nparts","jet_area","jet_angularity","jet_track_pt_0",
                        "jet_track_pt_1","jet_track_pt_2","jet_track_pt_3","jet_track_pt_4",
                        "jet_track_pt_5","jet_track_pt_6","jet_track_pt_7"]
    target = "jet_pt_truth" 
    all_features = features + [target]
    pdf = pdf[all_features].sample(frac=0.5, random_state=1).copy()
    return pdf
       
if __name__ == '__main__':
    # create argument parser
    parser = argparse.ArgumentParser(description='Train a model with a given set of hyperparameters')
    parser.add_argument('--hidden_nodes_start', type=int, nargs=3, help='number of nodes in each hidden layer', required=True)
    parser.add_argument('--hidden_nodes_end', type=int, nargs=3, help='number of nodes in each hidden layer', required=True)
    parser.add_argument('--output', type=str, help='output file name', required=True)
    args = parser.parse_args()
    hidden_nodes_start = args.hidden_nodes_start
    hidden_nodes_end = args.hidden_nodes_end
    output = args.output
    
    echo = f"hidden_nodes_start: {hidden_nodes_start}, hidden_nodes_end: {hidden_nodes_end}, output: {output}"
    # get pdf
    pdf = GetPDF()
    mses = []
    params = []
    for i in range(hidden_nodes_start[0], hidden_nodes_end[0]+1):
        for j in range(hidden_nodes_start[1], hidden_nodes_end[1]+1):
            for k in range(hidden_nodes_start[2], hidden_nodes_end[2]+1):
                hidden_nodes = [i,j,k]
                mse, param = TrainModel(pdf, hidden_nodes)
                mses.append(mse)
                params.append(param)
                print(f"mse: {mse}, params: {param}")

    data = {"mse": mses, "params": params}
    with open(output, 'w') as outfile:
        json.dump(data, outfile)
    
    print("Done")
    exit(0)

    
    
    
    
    
    
    
    
 
        