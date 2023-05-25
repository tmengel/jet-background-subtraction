#!/usr/bin/env python3

import pandas as pd
import numpy as np

import tensorflow as tf
import tensorflow.keras as keras
from keras import layers


def PlotModel(model, pdf, features, cutoff, filename, features_tex, target_tex):
    '''
    Plot the model
    '''
    
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as colors

    params = {'axes.labelsize': 18,
                'axes.linewidth' : 1.5,
                'font.size': 18,
                'font.family': 'Times New Roman',
                'mathtext.fontset': 'stix',
                'legend.fontsize': 20,
                'xtick.labelsize': 18,
                'ytick.labelsize': 20,
                'text.usetex': False,
                'lines.linewidth': 1,
                'lines.linestyle': ' ',
                'lines.markersize' : 6,
                'lines.markeredgewidth' : 1,
                'xtick.major.size' : 5,
                'xtick.minor.size' : 3,
                'xtick.major.width' : 2,
                'xtick.minor.width' : 1,
                'xtick.direction' : 'in',
                'ytick.major.size' : 5,
                'ytick.minor.size' : 3,
                'ytick.major.width' : 2,
                'ytick.minor.width' : 1,
                'ytick.direction' : 'in',
                'xtick.minor.visible' : True,
                'ytick.minor.visible' : True,
                'savefig.transparent': True,
                'errorbar.capsize': 1.5,
                }
    plt.rcParams.update(params)
    
    num_layers = len(model.layers)
    print('Number of layers: ', num_layers)
    layers = []
    for i in range(num_layers):
            layers.append(model.get_layer(index=i))

    layerX = {}
    layerY = {}
    layerOutput = {}
    layerProbability = {}
    layerWeights = {}
    layerBias = {}
    layerNactive = {}

    ymin = 1
    ymax = 2000
    
    nsamples = pdf.shape[0]
    X = pdf[features].values
    InputX = np.ones(X.shape[1])
    InputY = np.linspace(ymin,ymax,X.shape[1])
    InputProbability = np.ones(X.shape[1])
    for i in range(num_layers):
            layerWeights[i] = layers[i].get_weights()[0]
            layerBias[i] = layers[i].get_weights()[1]
            layerX[i] = np.ones(layerWeights[i].shape[1])*(i+2)
            if i == num_layers-1:
                    layerY[i] = np.array([(ymax-ymin)/2])
            else:
                    layerY[i] = np.linspace(ymin,ymax,layerWeights[i].shape[1])
            if i == 0:
                    layerOutput[i] = layers[i](X)
            else:
                    layerOutput[i] = layers[i](layerOutput[i-1])
                    
            layerProbability[i] = np.sum(layerOutput[i]>cutoff, axis=0)/nsamples
            layerNactive[i] = np.count_nonzero(layerProbability[i])

    fig = plt.figure(figsize=(15,20))
    ax = fig.add_subplot(111)        
    
    for ilayer in range(num_layers):
            for i in range(layerWeights[ilayer].shape[0]):
                    for j in range(layerWeights[ilayer].shape[1]):
                                    if ilayer == 0:
                                            if layerProbability[ilayer][j] > cutoff:
                                                    ax.plot([InputX[i]+0.02, layerX[ilayer][j]-0.02], [InputY[i], layerY[ilayer][j]], linewidth=0.5,
                                                            c=cm.Blues(np.abs(layerWeights[ilayer][i,j])/np.max(np.abs(layerWeights[ilayer][:,j]))), linestyle='-', alpha=0.7)
                                    else:
                                            if layerProbability[ilayer][j] > cutoff:
                                                    ax.plot([layerX[ilayer-1][i]+0.02, layerX[ilayer][j]-0.02], [layerY[ilayer-1][i], layerY[ilayer][j]], linewidth=0.5,
                                                            c=cm.Blues(np.abs(layerWeights[ilayer][i,j])/np.max(np.abs(layerWeights[ilayer][:,j]))), linestyle='-', alpha=0.7)
                                    
    ax.scatter(InputX, InputY, s=30, c=InputProbability, cmap='Reds', alpha=1.0, vmin=0, vmax=1)
    for i in range(num_layers):
            ax.scatter(layerX[i][layerProbability[i] > cutoff], layerY[i][layerProbability[i] > cutoff], s=30, c=layerProbability[i][layerProbability[i]>cutoff], cmap='Reds', alpha=1.0, vmin=0, vmax=1)
            
    ax.set_xlim(0.8, num_layers+1.2)
    xticks = [i+1 for i in range(num_layers+1)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(['Input']+['Layer '+str(i+1) for i in range(0,num_layers-1)]+['Output'])
    ax.set_ylim(ymin-200,ymax+200)
    ax.set_yticks(InputY)
    ax.set_yticklabels(features_tex)
    
    ax.cbar = plt.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=1), cmap='Reds'), ax=ax, location='top', pad=0.1, aspect=50, anchor=(0.5,-0.7))
    ax.cbar.set_label('Probability of Neuron Activation', rotation=0, labelpad=10)
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim()[0],ax.get_ylim()[1])
    
    ax2.set_yticks(layerY[num_layers-1])
    ax2.set_yticklabels(target_tex)
    ax2.cbar = plt.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=1), cmap='Blues'), ax=ax2, location='bottom', pad=0.1, aspect=50, anchor=(0.5,1.5))
    ax2.cbar.set_label('Normalized Weights', rotation=0, labelpad=10)
        
    active_nodes = []
    for i in range(num_layers):
            active_nodes.append(layerNactive[i])
    fig.savefig(f'plots/{filename}.pdf', bbox_inches='tight')
    return active_nodes
  
def MinimizeArch(pdf, features, target, cutoff, features_tex, target_tex):
    
    layernodes = {}
    mse_array = []
    
    inputlen = len(features)
    layer1 = layers.Dense(100,input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
    layer2 = layers.Dense(100,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
    layer3 = layers.Dense(100,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
    layer4 = layers.Dense(1,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
 
    model = keras.Sequential([layer1,layer2,layer3,layer4])
    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
    mse_initial = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
    mse_array.append(mse_initial)
    layernodes['iteration0'] = [100,100,100,1]
    
    active = PlotModel(model, pdf, features, cutoff, 'iteration0', features_tex, target_tex)    
    print("Initial Network: MSE: {mse} | Active Nodes: {active_nodes}".format(mse=mse_initial, active_nodes=active))
    
    mse_iter = 0 
    mse_best_idx = 0
    mse_best = mse_initial
    iteration = 0
    
    while True:
        #check for 0 nodes
        for act in active:
            if act == 0:
                active[active.index(act)] = 1
                l2_reg = 0.001
        
        # change l2 reg
        n_nodes_total = np.sum(active)
        if n_nodes_total < 50:
                l2_reg = 0.001
        else:
                l2_reg = 0.01
        
        if n_nodes_total < 10:
            break
        
        # new layer list
        layers_list = []
        # input layer    
        layers_list.append(layers.Dense(active[0],input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(l2_reg)))
        # hidden layers
        for act in active[1:]:
            layers_list.append(layers.Dense(act,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(l2_reg)))
       
        #new model
        model = keras.Sequential(layers_list)
        model.summary()
        model.compile(loss='mse', optimizer='adam', metrics=['mse'])    
       
        # save layer nodes
        layernodes[f'iteration{iteration+1}'] = active
        
        # fit model
        model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
        
        # evaluate model
        mse_iter = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
        mse_array.append(mse_iter)
     
        # plot model
        active = PlotModel(model, pdf, features, cutoff, f'iteration{iteration+1}', features_tex, target_tex)
        
        # print results
        print("Iteration {iter}: MSE: {mse} | Active Nodes: {active_nodes}".format(iter=iteration+1, mse=mse_iter, active_nodes=active))
        
        # iterate
        iteration += 1
        
        if mse_iter < mse_best:
            mse_best = mse_iter
            mse_best_idx = iteration-1
            continue
        
        mse_percent_change = (mse_iter - mse_best)/mse_best
        if mse_percent_change > 0.1:
            break
  
    print("Best Network 3-layer: MSE: {mse} | Active Nodes: {active_nodes}".format(mse=mse_best, active_nodes=layernodes[f'iteration{mse_best_idx}']))
    
    # now drop a layer
    print("Dropping a layer")
    n_nodes_total = np.sum(active)
    # new layer list
    layer1 = layers.Dense((n_nodes_total)+1,input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
    layer2 = layers.Dense(n_nodes_total,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
    layer3 = layers.Dense(1,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
    
    #new model
    model = keras.Sequential([layer1,layer2,layer3])
    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
    
    # evaluate model
    mse_initial = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
    
    # print model
    active = PlotModel(model, pdf, features, cutoff, f'iteration{iteration}', features_tex, target_tex)
    
    # save results
    mse_array.append(mse_initial)
    layernodes[f'iteration{iteration}'] = active
    
    # loop through iterations
    mse_best = mse_initial
    while True:
        #check for 0 nodes
        for act in active:
            if act == 0:
                active[active.index(act)] = 1
                l2_reg = 0.001
        
        # change l2 reg
        n_nodes_total = np.sum(active)
        if n_nodes_total < 25:
                l2_reg = 0.001
        else:
                l2_reg = 0.01
                
        if n_nodes_total < 10:
            break
        
        # new layer list
        layers_list = []
        # input layer    
        layers_list.append(layers.Dense(active[0],input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(l2_reg)))
        # hidden layers
        for act in active[1:]:
            layers_list.append(layers.Dense(act,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(l2_reg)))
       
        #new model
        model = keras.Sequential(layers_list)
        model.summary()
        model.compile(loss='mse', optimizer='adam', metrics=['mse'])    
       
        # save layer nodes
        layernodes[f'iteration{iteration+1}'] = active
        
        # fit model
        model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
        
        # evaluate model
        mse_iter = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
        mse_array.append(mse_iter)
     
        # plot model
        active = PlotModel(model, pdf, features, cutoff, f'iteration{iteration+1}', features_tex, target_tex)
        
        # print results
        print("Iteration {iter}: MSE: {mse} | Active Nodes: {active_nodes}".format(iter=iteration+1, mse=mse_iter, active_nodes=active))
        
        # iterate
        iteration += 1
        
        if mse_iter < mse_best:
            mse_best = mse_iter
            mse_best_idx = iteration-1
            continue
        
        mse_percent_change = (mse_iter - mse_best)/mse_best
        if mse_percent_change > 0.1:
            break


    print("Best 2-layer Network: MSE: {mse} | Active Nodes: {active_nodes}".format(mse=mse_best, active_nodes=layernodes[f'iteration{mse_best_idx}']))
    
    # now drop a layer
    print("Dropping a layer")
    n_nodes_total = np.sum(active)
    
    # new layer list
    layer1 = layers.Dense(2*(n_nodes_total)+1,input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
    layer2 = layers.Dense(1,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.01))
    
    #new model
    model = keras.Sequential([layer1,layer2])
    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
    
    # evaluate model
    mse_initial = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
    # plot model
    active = PlotModel(model, pdf, features, cutoff, f'iteration{iteration}', features_tex, target_tex)
    
    # save results
    mse_array.append(mse_initial)
    layernodes[f'iteration{iteration}'] = active
    
    mse_best = mse_initial
    while True:
        #check for 0 nodes
        for act in active:
            if act == 0:
                active[active.index(act)] = 1
                l2_reg = 0.001
        
        # change l2 reg
        n_nodes_total = np.sum(active)
        if n_nodes_total < 15:
                l2_reg = 0.001
        else:
                l2_reg = 0.01
                
        if n_nodes_total < 5:
            break
        
        # new layer list
        layers_list = []
        # input layer    
        layers_list.append(layers.Dense(active[0],input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(l2_reg)))
        # hidden layers
        for act in active[1:]:
            layers_list.append(layers.Dense(act,activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(l2_reg)))
       
        #new model
        model = keras.Sequential(layers_list)
        model.summary()
        model.compile(loss='mse', optimizer='adam', metrics=['mse'])    
       
        # save layer nodes
        layernodes[f'iteration{iteration+1}'] = active
        
        # fit model
        model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
        
        # evaluate model
        mse_iter = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
        mse_array.append(mse_iter)
     
        # plot model
        active = PlotModel(model, pdf, features, cutoff, f'iteration{iteration+1}', features_tex, target_tex)
        
        # print results
        print("Iteration {iter}: MSE: {mse} | Active Nodes: {active_nodes}".format(iter=iteration+1, mse=mse_iter, active_nodes=active))
        
        # iterate
        iteration += 1
        
        if mse_iter < mse_best:
            mse_best = mse_iter
            mse_best_idx = iteration-1
            continue
        
        mse_percent_change = (mse_iter - mse_best)/mse_best
        if mse_percent_change > 0.1:
            break
        
    print("Best 1-layer Network: MSE: {mse} | Active Nodes: {active_nodes}".format(mse=mse_best, active_nodes=layernodes[f'iteration{mse_best_idx}'])) 
    
    
    return mse_array, layernodes

def DNN(pdf, features, target, cutoff, features_tex, target_tex):
    
    inputlen = len(features)
    model = keras.Sequential(
            [
                layers.Dense(100,input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(100, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(50, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(1, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
            ])

    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
    mse = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
    layernodes = [100,100,50,1]
    
    # plot model
    active = PlotModel(model, pdf, features, cutoff, 'DNN', features_tex, target_tex)
    
    return mse, layernodes
    
def KSNN(pdf, features, target, cutoff, features_tex, target_tex):
    
    inputlen = len(features)
    model = keras.Sequential(
            [
                layers.Dense(100,input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(100, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(100, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(100, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(1, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
          ])

    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
    mse = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
    layernodes = [100,100,50,1]
    
    # plot model
    active = PlotModel(model, pdf, features, cutoff, 'kitchensink', features_tex, target_tex)
    
    return mse, layernodes

def SNN(pdf, features, target, cutoff, features_tex, target_tex):
    
    inputlen = len(features)
    model = keras.Sequential(
            [
                layers.Dense(inputlen+1,input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform'),
                layers.Dense(1, activation='relu',kernel_initializer='he_uniform'),
            ])

    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    model.fit(pdf[features], pdf[target], epochs=50, batch_size=64, verbose=0, validation_split=0.2)
    mse = model.evaluate(pdf[features], pdf[target], verbose=0)[1]
    layernodes = [inputlen+1,1]
    
    # plot model
    active = PlotModel(model, pdf, features, cutoff, 'SNN', features_tex, target_tex)
    
    return mse, layernodes

def MakeModelParameterPlot():
    
    print("Making Model Parameter Plot")
    file = "/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/root-files/AuAu_R04/test/AuAu_R04_50_test_sample.h5"
    pdf = pd.read_hdf(file, key='df')

    # minimize architecture
    features=['jet_angularity', 'jet_area', 'jet_pt_raw',  'jet_nparts', 'jet_track_pt_0']
    target='jet_pt_truth'
    features_tex = [r'$\lambda_{\alpha=1}$', r'$A_{\mathrm{jet}}$', r'$p_{\mathrm{T,raw}}$', r'$N_{\mathrm{track}}$', r'$p_{\mathrm{T,track}}^{0}$']
    target_tex = [r'$p_{\mathrm{T,truth}}$']
    cutoff = 0.0
    
    print("Minimizing Architecture")
    mse, layernodes = MinimizeArch(pdf = pdf, features = features, target = target, features_tex = features_tex, target_tex = target_tex, cutoff = cutoff)
    
    # calculate trainable parameters
    model_params = {}
    for key, value in zip(layernodes.keys(), mse):
        mod_params = {}
        trainable_params = 0
        for i in range(len(layernodes[key])):
            if i == 0:
                trainable_params += layernodes[key][i]*(len(features)+1)
            else:
                trainable_params += layernodes[key][i]*(layernodes[key][i-1]+1)
        
        mod_params['parameters'] = trainable_params
        mod_params['mse'] = value
        model_params[key] = mod_params            
    
    print("DNN")
    # DNN
    dnn_features = ["jet_pt_raw","jet_nparts","jet_area","jet_angularity","jet_track_pt_0",
                        "jet_track_pt_1","jet_track_pt_2","jet_track_pt_3","jet_track_pt_4",
                        "jet_track_pt_5","jet_track_pt_6","jet_track_pt_7"]
    dnn_features_tex = [r'$p_{\mathrm{T,raw}}$', r'$N_{\mathrm{track}}$', r'$A_{\mathrm{jet}}$', r'$\lambda_{\alpha=1}$', r'$p_{\mathrm{T,track}}^{0}$',
                        r'$p_{\mathrm{T,track}}^{1}$', r'$p_{\mathrm{T,track}}^{2}$', r'$p_{\mathrm{T,track}}^{3}$', r'$p_{\mathrm{T,track}}^{4}$',
                        r'$p_{\mathrm{T,track}}^{5}$', r'$p_{\mathrm{T,track}}^{6}$', r'$p_{\mathrm{T,track}}^{7}$']
    
    mse_dnn, layernodes_dnn = DNN(pdf = pdf, features = dnn_features, target = target, features_tex = dnn_features_tex, target_tex = target_tex, cutoff = cutoff)
    
    # calculate trainable parameters
    dnn_params = {}
    trainable_params=0
    for i in range(len(layernodes_dnn)):
        if i == 0:
            trainable_params += layernodes_dnn[i]*(len(dnn_features)+1)
        else:
            trainable_params += layernodes_dnn[i]*(layernodes_dnn[i-1]+1)
    
    dnn_params['parameters'] = trainable_params
    dnn_params['mse'] = mse_dnn
    model_params['DNN'] = dnn_params  
    
    print("Kitchen Sink")
    # kitchen sink
    kitchen_sink = ['jet_angularity',
             'jet_area', 'jet_average_track_pt', 
             'jet_nparts',
             'jet_pt_raw', 
             'jet_track_deltaR_0', 'jet_track_deltaR_1',
             'jet_track_energy_fraction_0', 'jet_track_energy_fraction_1', 
             'jet_track_pt_0', 'jet_track_pt_1', 'jet_track_pt_2', 'jet_track_pt_3',
             'jet_track_pt_4', 'jet_track_pt_5', 'jet_track_pt_6', 'jet_track_pt_7', 'jet_track_pt_8', 'jet_track_pt_9',
             'jet_track_pt_kurtosis',  'jet_track_pt_skewness', 
             'jet_track_pt_variance', 'median_pt_over_area', 'median_pt_over_npart', 'random_cone_nparts', 'random_cone_pt']
    
    kitchen_sink_tex = [r'$\lambda_{\alpha=1}$',
                r'$A_{\mathrm{jet}}$', r'$\langle p_{\mathrm{T,track}} \rangle$',
                r'$N_{\mathrm{track}}$',
                r'$p_{\mathrm{T,raw}}$',
                r'$\Delta R_{0}$', r'$\Delta R_{1}$',
                r'$f_{E,0}$', r'$f_{E,1}$',
                r'$p_{\mathrm{T,track}}^{0}$', r'$p_{\mathrm{T,track}}^{1}$', r'$p_{\mathrm{T,track}}^{2}$', r'$p_{\mathrm{T,track}}^{3}$',
                r'$p_{\mathrm{T,track}}^{4}$', r'$p_{\mathrm{T,track}}^{5}$', r'$p_{\mathrm{T,track}}^{6}$', r'$p_{\mathrm{T,track}}^{7}$', r'$p_{\mathrm{T,track}}^{8}$', r'$p_{\mathrm{T,track}}^{9}$',
                r'$\kappa_{p_{\mathrm{T,track}}}$', r'$\gamma_{p_{\mathrm{T,track}}}$',
                r'$\sigma_{p_{\mathrm{T,track}}}$', r'$\langle p_{\mathrm{T,track}}/A_{\mathrm{jet}} \rangle$', r'$\langle p_{\mathrm{T,track}}/N_{\mathrm{track}} \rangle$', r'$N_{\mathrm{track}}^{\mathrm{cone}}$', r'$p_{\mathrm{T,track}}^{\mathrm{cone}}$']
    
    mse_kitchensink, layernodes_kitchensink = KSNN(pdf = pdf, features = kitchen_sink, target = target, features_tex = kitchen_sink_tex, target_tex = target_tex, cutoff = cutoff)
    
    # calculate trainable parameters
    ks_params = {}
    trainable_params=0
    for i in range(len(layernodes_kitchensink)):
        if i == 0:
            trainable_params += layernodes_kitchensink[i]*(len(kitchen_sink)+1)
        else:
            trainable_params += layernodes_kitchensink[i]*(layernodes_kitchensink[i-1]+1)
    
    ks_params['parameters'] = trainable_params
    ks_params['mse'] = mse_kitchensink
    model_params['KitchenSink'] = ks_params
    
    print("SNN")
    # SNN 
    mse_snn, layernodes_snn = SNN(pdf = pdf, features = features, target = target, features_tex = features_tex, target_tex = target_tex, cutoff = cutoff)
    
    # calculate trainable parameters
    snn_params = {}
    trainable_params=0
    for i in range(len(layernodes_snn)):
        if i == 0:
            trainable_params += layernodes_snn[i]*(len(features)+1)
        else:
            trainable_params += layernodes_snn[i]*(layernodes_snn[i-1]+1)
    
    snn_params['parameters'] = trainable_params
    snn_params['mse'] = mse_snn
    model_params['SNN'] = snn_params
    
    import matplotlib.pyplot as plt
    params = {'axes.labelsize': 18,
                'axes.linewidth' : 1.5,
                'font.size': 18,
                'font.family': 'Times New Roman',
                'mathtext.fontset': 'stix',
                'legend.fontsize': 20,
                'xtick.labelsize': 18,
                'ytick.labelsize': 20,
                'text.usetex': False,
                'lines.linewidth': 1,
                'lines.linestyle': ' ',
                'lines.markersize' : 6,
                'lines.markeredgewidth' : 1,
                'xtick.major.size' : 5,
                'xtick.minor.size' : 3,
                'xtick.major.width' : 2,
                'xtick.minor.width' : 1,
                'xtick.direction' : 'in',
                'ytick.major.size' : 5,
                'ytick.minor.size' : 3,
                'ytick.major.width' : 2,
                'ytick.minor.width' : 1,
                'ytick.direction' : 'in',
                'xtick.minor.visible' : True,
                'ytick.minor.visible' : True,
                'savefig.transparent': True,
                'errorbar.capsize': 1.5,
                }
    plt.rcParams.update(params)
    # plot results
    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_ylabel(r'$\mathrm{Mean \ Squared \ Error} (GeV^2)$')
    ax.set_xlabel(r'$\mathrm{Number \ of \ Trainable \ Parameters}$')
    ax.set_title(r'$\mathrm{Model \ Performance \ Curve}$')
    ax.set_xscale('log')
    
    parameter_array = []
    mse_array = []
    for key in model_params:
        parameter_array.append(model_params[key]['parameters'])
        mse_array.append(model_params[key]['mse'])
        
    ax.plot(parameter_array, mse_array, 'o', label = r'$\mathrm{Model \ Performance}$', color = 'black')
    
    # highlight SNN, DNN, and Kitchen Sink
    ax.plot(model_params['SNN']['parameters'], model_params['SNN']['mse'], 'o', label = r'$\mathrm{Shallow Neural Network}$', color = 'red')
    ax.plot(model_params['DNN']['parameters'], model_params['DNN']['mse'], 'o', label = r'$\mathrm{Deep Neural Network}$', color = 'blue')
    ax.plot(model_params['KitchenSink']['parameters'], model_params['KitchenSink']['mse'], 'o', label = r'$\mathrm{Kitchen \ Sink}$', color = 'green')
    
    ax.legend(loc = 'best', frameon = False)
    fig.savefig('plots/model_optimization.pdf', bbox_inches='tight')
    
if __name__ == '__main__':
    print("Starting")    
    MakeModelParameterPlot()
    print("Finished")
   
    
    
    
    
    
    
    
 
        