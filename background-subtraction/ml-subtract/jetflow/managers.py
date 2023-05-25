from __future__ import absolute_import

import os
import json 

from jetflow.base import *
from jetflow.utils import *
    
__all__ = ['DataManager', 'ModelManager']

## Data Manager class

class DataManager(JFBase):
    
    '''
        Data Manager class handles loading and preprocessing of data
        from ROOT files. Will also handle splitting data into training
        and testing sets. All data stored in global data directory after 
        it has been loaded and preprocessed for the first time.
            
    '''
    def __init__(self):
        super(DataManager,self).__init__()
        self.ConfigPath = os.path.join(self.DataDirectory,'DataConfig.json') 
        self.Config = None
        self.Datasets = []
        self.DatasetObjects = {}
        
        # Get the path to the data configuration file
        if os.path.isfile(self.ConfigPath):
            # If the file exists, load it
            conf = {}
            with open(self.ConfigPath) as f:
                conf = json.load(f)
                
            self.Config = conf
            self.Datasets = self.Config['datasets']
            
        else:
            self.Config = {}
            self.Config['datasets'] = []
            self.Datasets = self.Config['datasets']
        
        self.SaveConfig()
        
    def UpdateConfig(self):
        for ds in self.Datasets:
            if ds not in self.Config['datasets']:
                self.Config['datasets'].append(ds)
    
        self.SaveConfig()
    
    def SaveConfig(self):
        with open(self.ConfigPath,'w') as f:
            json.dump(self.Config,f,indent=4)
        
    def Add(self, name, filename, overwrite=False):
        '''
            Add a dataset to the data manager. This will create a directory for the dataset
            in the data directory and store the configuration file there.
        '''
        if name in self.Datasets:
            if not overwrite:
                print('Dataset already exists. Use overwrite=True to overwrite.')
                return
            else:
                self.Datasets.remove(name)
                
        DS  = Dataset(name, filename, overwrite=overwrite)
        self.Datasets.append(name)
        self.DatasetObjects.update({name:DS})
        self.UpdateConfig()
      
    def Get(self, name):
        if name not in self.Datasets:
            raise ValueError('Dataset not found: ' + name)
        return self.DatasetObjects[name]
        
    def Describe(self):
        '''
            List the datasets that are available.
        '''
        for ds in self.Datasets:
            print(ds)
    
    def PrintConfig(self):
        '''
            Print the configuration file.
        '''
        print(json.dumps(self.Config,indent=4))
    
    def Remove(self, name):
        '''
            Remove a dataset from the data manager.
        '''
        if name not in self.Datasets:
            raise ValueError('Dataset not found: ' + name)
        self.Datasets.remove(name)
        self.DatasetObjects['datasets'].pop(name)
        self.UpdateConfig()
    
    def __call__(self, name):
        self.UpdateConfig()
        return self.Get(name)
        
## Model Manager class

class ModelManager(JFBase):
    
    def __init__(self):
        super(ModelManager,self).__init__()
        self.ConfigPath = os.path.join(self.ModelDirectory,'ModelConfig.json') 
        self.Config = None
        self.Models = []
        self.ModelObjects = {}
        
         # Get the path to the data configuration file
        self.ConfigPath = os.path.join(self.ModelDirectory,'ModelConfig.json')
        if os.path.isfile(self.ConfigPath):
            
            conf = {}
            with open(self.ConfigPath) as f:
                conf = json.load(f)
            
            self.Config = conf
            self.Models = self.Config['models']
        
        else:
            self.Config = {}
            self.Config['models'] = []
            self.Models = self.Config['models']
        
        self.SaveConfig()
     
    def UpdateConfig(self):
        self.Config.update({'models':self.Models})
        self.SaveConfig()
        
    def SaveConfig(self):
        with open(self.ConfigPath,'w') as f:
            json.dump(self.Config,f,indent=4)
        
    def Add(self, name, type, features, target, retrain=False):
        '''
            Add a new model to the model manager. This will create a new model configuration
        '''
        if name in self.Models and not retrain:
            print('Model already exists. Set overwrite=True to overwrite.')
            return
        elif name in self.Models and retrain:
            print('Overwriting model.')
            self.Models.remove(name)
        
        Mod = Model(name, type, features, target, retrain)
        self.Models.append(name)
        self.ModelObjects.update({name:Mod})
        self.UpdateConfig()
        
    def GetModel(self, name):
        if name in self.Models:
            return self.ModelObjects[name]
        else:
            print('Model not found.')
            return None
    
    def Describe(self):
        '''
            List the datasets that are available.
        '''
        for ds in self.Models:
            print(ds)
    
    def PrintConfig(self):
        '''
            Print the configuration file.
        '''
        print(json.dumps(self.Config,indent=4))
            
        

