# Author: Tanner Mengel
# JetFlow is a class that contains dataset and model management for machine learning
# based background subtraction. 

from __future__ import absolute_import

import os
import json 

from jetflow.base import *
from jetflow.utils import *
from jetflow.managers import *
            
            
# __all__ = ['JetFlow']

class JF(JFBase):
    '''
        JetFlow class handles dataset and model management for machine learning
        based background subtraction. 
    '''
    
    def __init__(self, InputConfig):
        super(JF,self).__init__()
        self.DataManager = DataManager()
        self.ModelManager = ModelManager()
        self.Config = None
        
        # make sure config file exists
        if not os.path.isfile(InputConfig):
            raise IOError('Input config file does not exist')
        # process the input config
        # nessisary to make sure that the config is in the correct format
        # need dataset name, source files, model types
        conf = {}
        with open(InputConfig) as f:
            conf = json.load(f)
        keys = conf.keys()
        if 'Datasets' not in keys:
            raise ValueError('Input config file does not contain Datasets key')
        if 'Models' not in keys:
            raise ValueError('Input config file does not contain Models key')
        
        # check that all nessisary keys are present
        dataconf = conf['Datasets']
        dat_nessisary_keys = ['Name', 'Path', 'Overwrite', 'Split']
        for k in dat_nessisary_keys:
                if k not in dataconf.keys():
                    raise ValueError('Input config file does not contain nessisary keys: {}'.format(k))
       
        modelconf = conf['Models']
        mod_nessisary_keys = ['Name', 'Type', 'Features', 'Target', 'Retrain', 'Epochs', 'Batch Size', 'Optimizer', 'Loss']
        for mod in modelconf:
            modconf_keys = modelconf[mod].keys()
            for k in mod_nessisary_keys:
                if k not in modconf_keys:
                    raise ValueError('Input config file does not contain nessisary keys: {}'.format(k))
    
        # if all nessisary keys are present, then we can load the config
        self.Config = conf
        
        
        # add datasets
       
        self.DataManager.Add(self.Config['Datasets']['Name'], self.Config['Datasets']['Path'], self.Config['Datasets']['Overwrite'])
        self.DataManager.Get(self.Config['Datasets']['Name']).Split(self.Config['Datasets']['Split'])
        
        # add models
        for mod in self.Config['Models']:
            self.ModelManager.Add(self.Config['Models'][mod]['Name'],
                                  self.Config['Models'][mod]['Type'],
                                  self.Config['Models'][mod]['Features'],
                                  self.Config['Models'][mod]['Target'],
                                  self.Config['Models'][mod]['Retrain'])
            
            self.ModelManager.GetModel(self.Config['Models'][mod]['Name']).Compile(self.Config['Models'][mod]['Optimizer'], self.Config['Models'][mod]['Loss'])
        
        print('JetFlow initialized')
        
    def Train(self):
        '''
            Train all models
        '''
     
        dataset = self.DataManager.Get(self.Config['Datasets']['Name'])
        train = dataset.Load('train')
        for mod in self.Config['Models']:
            print('Training model: {}'.format(self.Config['Models'][mod]['Name']))
            model = self.ModelManager.GetModel(self.Config['Models'][mod]['Name'])
            x = train[model.features]
            y = train[model.target]
            model.Train(x, y, epochs=self.Config['Models'][mod]['Epochs'], batch_size=self.Config['Models'][mod]['Batch Size'], val_split=0.2)
            model.Save()
        
        self.SaveModelsToTFile()
    
    def Eval(self):
        '''
            Evaluate all models
        '''
        dataset = self.DataManager.Get(self.Config['Datasets']['Name'])
        test = dataset.Load('test')
        for mod in self.Config['Models']:
            print('Evaluating model: {}'.format(self.Config['Models'][mod]['Name']))
            model = self.ModelManager.GetModel(self.Config['Models'][mod]['Name'])
            x = test[model.features]
            test[f'jet_pt_{model.type}_corrected'] = model.model.predict(x, verbose=0)[:,0]
            
        print('Saving results')
        # save the results
        output_root_file_name = os.path.join(self.ResultsDirectory, f'{self.Config["Datasets"]["Name"]}_results.root')
        WriteTFile(test, output_root_file_name)
        print('Saving results to hdf5')
        # save hdf5 file
        test.to_hdf(os.path.join(self.ResultsDirectory, f'{self.Config["Datasets"]["Name"]}_results.h5'), key='df', mode='w')
        print('Done')
            
    def SaveModelsToTFile(self):
        '''
            Save all models to a root file
        '''
        import ROOT
        import numpy as np
        for mod in self.Config['Models']:
            
            model = self.ModelManager.GetModel(self.Config['Models'][mod]['Name'])
            filename = model.path.replace('.h5', '.root')
            
            # get the model weights and biases and activation functions           
            weights = model.model.get_weights()
            num_layers = len(weights)//2
            weight_array = []
            bias_array = []
            activation_functions = []
            for i in range(num_layers):
                activation_functions.append(model.model.get_layer(index=i).activation.__name__)
                weight_array.append(weights[2*i])
                bias_array.append(weights[2*i+1])
            
            # create the root file
            fout = ROOT.TFile(filename, 'RECREATE')
            tree_architecture = ROOT.TTree("arch","arch")
            tree_activations = ROOT.TTree("activations","activations")
            tree_inputs = ROOT.TTree("inputs","inputs")
            tree_targets = ROOT.TTree("targets","targets")
            
            # achitecture tree variables
            n_activations = np.zeros(1,dtype=np.int32)
            n_inputs = np.zeros(1,dtype=np.int32)
            n_targets = np.zeros(1,dtype=np.int32)
            model_type = np.zeros(1,dtype='|S20')

            # architecture tree
            tree_architecture.Branch("n_layers",n_activations,"n_layers/I")
            tree_architecture.Branch("n_inputs",n_inputs,"n_inputs/I")
            tree_architecture.Branch("n_targets",n_targets,"n_targets/I")
            tree_architecture.Branch("model_type",model_type,"model_type/C")
            model_type[0] = model.type
            n_activations[0] = num_layers
            n_inputs[0] = len(model.features)
            n_targets[0] = len(model.target)
            tree_architecture.Fill()
            
            # activations tree variables
            acts = np.zeros(1,dtype='|S20')
            # activations tree
            tree_activations.Branch("activations",acts,"activations/C")
            for i in range(len(activation_functions)):
                acts[0] = activation_functions[i]
                tree_activations.Fill()
                
            # input tree variables
            inputs = np.zeros(1,dtype='|S20')
            input_types = np.zeros(1,dtype='|S20')
            # input tree
            tree_inputs.Branch("inputs",inputs,"inputs/C")
            tree_inputs.Branch("input_types",input_types,"input_types/C")
            for i in range(len(model.features)):
                inputs[0] = model.features[i]
                if model.features[i] == "jet_nparts":
                    input_types[0] = 'Int_t'
                else:
                    input_types[0] = 'Float_t'
                tree_inputs.Fill()
                
            # target tree variables
            targets = np.zeros(1,dtype='|S20')
            target_types = np.zeros(1,dtype='|S20')
            # target tree
            tree_targets.Branch("targets",targets,"targets/C")
            tree_targets.Branch("target_types",target_types,"target_types/C")
            for i in range(len(model.target)):
                targets[0] = model.target[i]
                target_types[0] = 'Float_t'
                tree_targets.Fill()

            # weights and biases
            for i in range(len(weight_array)):
                th2d = ROOT.TH2D(f'weight{i}',f'weight{i}',weight_array[i].shape[0],0,weight_array[i].shape[0],weight_array[i].shape[1],0,weight_array[i].shape[1])
                th1d = ROOT.TH1D(f'bias{i}',f'bias{i}',bias_array[i].shape[0],0,bias_array[i].shape[0])
                for j in range(weight_array[i].shape[0]):
                    for k in range(weight_array[i].shape[1]):
                        th2d.SetBinContent(j+1,k+1,weight_array[i][j][k])
                for j in range(bias_array[i].shape[0]):
                    th1d.SetBinContent(j+1,bias_array[i][j])
                th2d.Write()
                th1d.Write()
                
                
            fout.Write()
            fout.Close()
                    
           
                
        
