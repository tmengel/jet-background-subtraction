# Base directory for the jetflow package
from __future__ import absolute_import, division, print_function

import os
import warnings
import pandas as pd
import json 
from jetflow.utils import *

try:
    from sklearn.model_selection import train_test_split
except ImportError as e:
    print('Could not import - ' + str(e))

###########################################################
## Global Properties of the jetflow package
###########################################################

__all__ = ['JFBase', 'Dataset', 'Model']

class JFBase(object):
    '''
        Global variables for the jetflow package.
    '''
    def __init__(self):
        self.data_dir = self.GetDataDirectory()
        self.model_dir = self.GetModelDirectory()
        self.results_dir = self.GetResultsDirectory()
        
    @property
    def DataDirectory(self):
        return self.data_dir if self.data_dir else self.GetDataDirectory()
    
    @property
    def ModelDirectory(self):
        return self.model_dir if self.model_dir else self.GetModelDirectory()
    
    @property
    def ResultsDirectory(self):
        return self.results_dir if self.results_dir else self.GetResultsDirectory()
    
    @staticmethod
    def GetDataDirectory():
        '''
            Find the data directory. This is the directory where all data is stored.
            If data directory does not exist, create it. The data directory should be
            in the same directory as the file that imports the jetflow package.
        '''
        # Get the directory of file that imports the jetflow package
        dir_path = os.path.dirname(os.path.realpath(__file__))
        # Get the parent directory of this file
        parent_dir = os.path.dirname(dir_path)
        # Get the data directory
        data_dir = os.path.join(parent_dir,'data')
        # Check if the data directory exists
        if os.path.isdir(data_dir):
            return data_dir
        else:
            # If the data directory does not exist, create it
            os.mkdir(data_dir)
            return data_dir

    @staticmethod
    def GetModelDirectory():
        '''
            Find the model directory. This is the directory where all models are stored.
            If model directory does not exist, create it. The model directory should be
            in the same directory as the file that imports the jetflow package.
        '''
        # Get the directory of file that imports the jetflow package
        dir_path = os.path.dirname(os.path.realpath(__file__))
        # Get the parent directory of this file
        parent_dir = os.path.dirname(dir_path)
        # Get the data directory
        data_dir = os.path.join(parent_dir,'models')
        # Check if the data directory exists
        if os.path.isdir(data_dir):
            return data_dir
        else:
            # If the data directory does not exist, create it
            os.mkdir(data_dir)
            return data_dir
        
    @staticmethod
    def GetResultsDirectory():
        '''
            Find the results directory. This is the directory where all results are stored.
            If results directory does not exist, create it. The results directory should be
            in the same directory as the file that imports the jetflow package.
        '''
        # Get the directory of file that imports the jetflow package
        dir_path = os.path.dirname(os.path.realpath(__file__))
        # Get the parent directory of this file
        parent_dir = os.path.dirname(dir_path)
        # Get the data directory
        data_dir = os.path.join(parent_dir,'results')
        # Check if the data directory exists
        if os.path.isdir(data_dir):
            return data_dir
        else:
            # If the data directory does not exist, create it
            os.mkdir(data_dir)
            return data_dir

class Dataset(JFBase):
    '''
        Dataset class stores information about a dataset and methods for loading
    '''
    def __init__(self, name, filename, overwrite=False):
        super(Dataset,self).__init__()
        self.name = name
        self.src = filename
        self.overwrite = overwrite
        self.train = None
        self.test = None
        self.features = None
        self.types = None
        self.tree = None
        self.nentries = None
        self.ntrain = None
        self.ntest = None
        self.path = None
        self.config_file = name+'.json'
        self.config = {}
        
        # check if filename is a existing file
        if not os.path.isfile(filename):
            raise ValueError('Src file is not a file: ' + filename)

        self.Configure()
        
        
    def Configure(self):
        # check if the dataset has already been loaded
        if self.name in os.listdir(self.DataDirectory):
            self.path = os.path.join(self.DataDirectory, self.name)
            
            # check if the dataset has configuration file
            if self.config_file in os.listdir(self.path):
                # check if overwrite is True
                if self.overwrite:
                    # if so, reconfigure the dataset
                    self.Create()
                elif not self.overwrite:
                    # if not, load the configuration
                    self.LoadConfig()
                
        # if not, create a directory for the dataset
        else:
            self.Create()
        
    def Create(self):
        self.path = os.path.join(self.DataDirectory, self.name)
        if not os.path.isdir(self.path):
            os.mkdir(self.path)
        if not os.path.isfile(os.path.join(self.path, self.config_file)):
            os.mknod(os.path.join(self.path, self.config_file))
        
        file_config = GetTFileConfig(self.src)
        self.tree = file_config['tree']
        self.nentries = file_config['entries']
        self.features = file_config['features']
        self.types = file_config['types']
        self.Update()
        
    def Update(self):
        self.config.update({'name':self.name})
        self.config.update({'src':self.src})
        self.config.update({'tree':self.tree})
        self.config.update({'path':self.path})
        self.config.update({'nentries':self.nentries})
        self.config.update({'features':self.features})
        self.config.update({'types':self.types})
        self.config.update({'ntrain':self.ntrain})
        self.config.update({'ntest':self.ntest})
        self.config.update({'train':self.train})
        self.config.update({'test':self.test})
        self.SaveConfig()
    
    def SaveConfig(self):
        # Save the configuration file in the path directory
        with open(os.path.join(self.path, self.config_file),'w') as f:
            json.dump(self.config, f, indent=4)
    
    def LoadConfig(self):
        # Load the configuration file from the path directory
        with open(os.path.join(self.path, self.config_file),'r') as f:
            self.config = json.load(f)
       
    def Load(self, type=None, features=None):
        
        if features is None:
            features = self.features
        for feature in features:
            if feature not in self.features:
                raise ValueError('Feature not in dataset: ' + feature)
            
        if type == 'train':
            if self.train is None:
                raise ValueError("Dataset has not been split yet. Use 'split' method first.")
            else:
                return LoadTFile(self.train, self.tree, features)
            
        elif type == 'test':
            if self.test is None:
                raise ValueError("Dataset has not been split yet. Use 'split' method first.")
            else:
                return LoadTFile(self.test, self.tree, features)
            
        elif type is None:
            return LoadTFile(self.src, self.tree, features)
        
        else:
            raise ValueError('Invalid type: ' + type)
        
    def Split(self, test_size=0.2, random_state=42):
        # load source file
        pdf = LoadTFile(self.src, self.tree, self.features)
        
        # split the dataset into train and test
        train_pdf, test_pdf = train_test_split(pdf, test_size=test_size, random_state=random_state)
       
        # remove any rows with NaN values
        train_pdf.dropna(inplace=True)
        test_pdf.dropna(inplace=True)
        
        # get the number of entries in the train and test datasets
        self.train = os.path.join(self.path, 'train.root')
        self.test = os.path.join(self.path, 'test.root')
        self.ntrain = len(train_pdf)
        self.ntest = len(test_pdf)
        
        # convert the pandas dataframes to root files and save them
        WriteTFile(train_pdf, self.train, self.tree)
        WriteTFile(test_pdf, self.test, self.tree)
        
        # update the configuration file
        self.Update()
        
    def Desribe(self):
        for key in self.config.keys():
            print(key + ' : ' + str(self.config[key]))
    
    def PrintConfig(self):
        return self.config
    
    def __call__(self, *args):
        self.UpdateConfig()
        return self.Load(*args)   
    
class Model(JFBase):
    '''
        Model class stores information about a models and methods for loading
    '''
    
    def __init__(self, name, type, features, target, retrain=False):
        super(Model,self).__init__()
        self.name = name
        self.type = type
        self.features = features
        self.target = target
        self.retrain = retrain
        self.model = None
        self.trained = False
        self.history = None
        self.loss = None
        self.optimizer = None
        self.epoch = None
        self.batch_size = None
        self.path = os.path.join(self.ModelDirectory, f'{self.name}_{self.target}.h5') 
        self.config_file = os.path.join(self.ModelDirectory, f'{self.name}_{self.target}.json')
        self.config = {}
        self.Configure()
                  
    def Configure(self):
        if os.path.isfile(self.path) and os.path.isfile(self.config_file):
            if self.retrain:
                warnings.warn('Model already exists. Retraining model.')
                self.Create()
            else:
                self.LoadConfig()
        else:
            self.Create()
                   
    def Create(self):
        
        # create the model file if it does not exist
        if not os.path.isfile(self.path):
            os.mknod(self.path)
        # create the configuration file if it does not exist
        if not os.path.isfile(self.config_file):
            os.mknod(self.config_file)
        
        # get the model
        self.model = GetModel(self.type, self.features)
        
        # update the configuration file
        self.Update()
        
    def Update(self):
        # update the configuration file
        self.config.update({'name':self.name})
        self.config.update({'type':self.type})
        self.config.update({'features':self.features})
        self.config.update({'target':self.target})
        self.config.update({'path':self.path})
        self.config.update({'trained':self.trained})
        self.config.update({'loss':self.loss})
        self.config.update({'optimizer':self.optimizer})
        self.config.update({'epoch':self.epoch})
        self.config.update({'batch_size':self.batch_size})
        # save the configuration file
        self.SaveConfig()
    
    def SaveConfig(self):
        # Save the configuration file in the path directory
        with open(self.config_file,'w') as f:
            json.dump(self.config, f, indent=4)
    
    def LoadConfig(self):
        # Load the configuration file from the path directory
        with open(os.path.join(self.path, self.config_file),'r') as f:
            self.config = json.load(f)
        self.name = self.config['name']
        self.type = self.config['type']
        self.features = self.config['features']
        self.target = self.config['target']
        self.path = self.config['path']
        self.trained = self.config['trained']
        self.loss = self.config['loss']
        self.optimizer = self.config['optimizer']
        self.epoch = self.config['epoch']
        self.batch_size = self.config['batch_size']
        self.model = LoadModel(self.path)
        self.Summary()

    def Desribe(self):
        for key in self.config.keys():
            print(key + ' : ' + str(self.config[key]))
        return self.config
    
    def PrintConfig(self):
        return self.config

    def Compile(self, optimizer='adam', loss='mse', metrics=['mse']):
        self.optimizer = optimizer
        self.loss = loss
        self.metrics = metrics
        self.Update()
        self.model.compile(optimizer=optimizer, loss=loss, metrics=metrics)
    
    def Summary(self):
        self.model.summary()
    
    def Train(self, x, y, val_split=0.2, epochs=100, batch_size=32, callbacks=None):
        self.history = self.model.fit(x, y, validation_split=val_split, epochs=epochs, batch_size=batch_size, callbacks=callbacks, verbose=0)
        self.trained = True
        self.epoch = epochs
        self.batch_size = batch_size
        self.Update()
        return self.history
    
    def Save(self):
        SaveModel(self.model, self.path)
    
    def Predict(self, x):
        return self.model.predict(x)
    
    def __call__(self):
        self.Update()
