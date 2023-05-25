#!/usr/bin/env python3
import pandas as pd
import numpy as np

import tensorflow as tf
import tensorflow.keras as keras
from keras import layers

import argparse
import warnings
    
def GetPrior(filename):
    '''
    Load prior file
    '''
    import json
    with open(filename) as json_file:
        data = json.load(json_file)
    return data  

def CopyDataFile(jobID,filename):
    ''''
    Copy data file from EOS to local directory
    '''
    
    import os
    import shutil
    
    # check if file exists
    if not os.path.isfile(filename):
        warnings.warn(f"File {filename} does not exist")
        exit(1)
    
    new_filename = f"tmpfiles/{jobID}.h5"
    # check if file already exists
    if os.path.isfile(new_filename):
        # if it does, delete it
        os.remove(new_filename)
    
    # copy file
    shutil.copyfile(filename, new_filename)
    
    return new_filename

def CopyPriorFile(jobID,filename):
    '''
    Copy prior file from EOS to local directory
    '''
    import os
    import shutil
       
    # check if file exists
    if not os.path.isfile(filename):
        warnings.warn(f"File {filename} does not exist")
        exit(1)
    
    new_filename = f"tmpfiles/{jobID}.json"
    
    # check if file already exists
    if os.path.isfile(new_filename):
        # if it does, delete it
        os.remove(new_filename)
    
    # copy file
    shutil.copyfile(filename, new_filename)
    
    return new_filename
    
def SamplePrior(prior):
   
    import time 
    currtime = time.time()
    seed = int((currtime%int(currtime))*100000)
    np.random.seed(seed)
    
    # sample prior
    selected_runs = {}
    for i in range(10):
        runs = {}
        # nvariables = [np.random.choice([i for i in range(1,prior['nVariables']+1)],1, p=prior['nVariables_weight'])[0]]
        nvariables = prior['nVariables']
        runs['variables_selected'] =[ "jet_angularity", "jet_area", "jet_nparts",  "jet_pt_raw",  "jet_track_pt_0"]
        nlayers = np.random.choice([i for i in range(1,prior['MaxLayers']+1)],1, p=prior['nLayers'])[0]
        nodes = []
        for j in range(nlayers):
            nodes.append(np.random.choice([i for i in range(1,prior['MaxNodes']+1)],1, p=prior[f'layer{j+1}'])[0]*1.0)
        runs['nodes'] = nodes
        selected_runs[i] = runs
    
    return selected_runs

def Run(run, pdf, target):
    '''
    reads selected parameters, creates model, evaluates model and returns the results
    '''
    variables = run['variables_selected']
    nodes = run['nodes']
    nlayers = len(nodes)
    inputdim = len(variables)
    
    layer_list = []
    layer_list.append(layers.Dense(int(nodes[0]), activation='relu', input_shape=(inputdim,),kernel_initializer='he_uniform'))
    layer_list.append(layers.Dropout(0.1))
    for i in range(1,nlayers):
        layer_list.append(layers.Dense(int(nodes[i]), activation='relu',kernel_initializer='he_uniform'))
        layer_list.append(layers.Dropout(0.1))
    layer_list.append(layers.Dense(1, activation='relu',kernel_initializer='he_uniform'))
    
    model = keras.Sequential(layer_list)
    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    model.fit(pdf[variables], pdf[target], epochs=25, batch_size=100, verbose=0, validation_split=0.2)
    loss, mse = model.evaluate(pdf[variables], pdf[target], verbose=0)
    
    trainable_parms = 0
    trainable_parms += (inputdim+1)*int(nodes[0])
    for i in range(1,nlayers):
        trainable_parms += (int(nodes[i-1])+1)*int(nodes[i])
    trainable_parms += (int(nodes[nlayers-1])+1)*1
    
    results = {}
    results['loss'] = [mse]
    results['trainable_parms'] = [trainable_parms]
    results['variables'] = [variables]
    results['nodes'] = [nodes]
    
    return results

def SaveResults(results, jobID):
    '''
    save results to file
    '''
    import json
    with open(f"runs/{jobID}.json", 'w') as outfile:
        json.dump(results, outfile, indent=4)

def DeleteTmpFiles(jobID):
    '''
    delete temporary files
    '''
    import os
    os.remove(f"tmpfiles/{jobID}.h5")
    os.remove(f"tmpfiles/{jobID}.json")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Hyperparameter Sampling Module')
    parser.add_argument('-j', '--job', required=True, help='job ID', type=str)
    args = parser.parse_args()
    
    src_filename= "/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/root-files/AuAu_R04/test/AuAu_R04_50_test_sample.h5"
    src_prior = "prior.json"

    # get prior
    prior = GetPrior(CopyPriorFile(args.job,src_prior))  
    # open hdf5 file
    pdf = pd.read_hdf(CopyDataFile(args.job,src_filename), key='df')
    
    # sample prior
    selected_runs = SamplePrior(prior)
    
    results = {}
    for run in selected_runs:
        print(f"Run {run}")
        results[run] = Run(selected_runs[run], pdf, prior['target'])
        # save results
        SaveResults(results, args.job)
    
    # save final results
    SaveResults(results, args.job)
    
    # delete copied files
    DeleteTmpFiles(args.job)
    
    print("Done")
    
    
    
    
    
    
    
    
    
    
   
    
    
    
    
    
    
    
 
        