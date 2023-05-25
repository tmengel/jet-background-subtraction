# root utils for loading and writing ROOT files
# reading trees and histograms

from __future__ import absolute_import, division, print_function

from jetflow.models import *

try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    
except ImportError as e:
    print('Could not import - ' + str(e))

__all__ = ['GetModel','SaveModel','LoadModel']


def GetModel(type, features):
    
    if type == "dnn":
        return dnn(features)
    elif type == "snn":
        return snn(features)   
    else:
        raise ValueError('Model type not recognized: ' + type)

def SaveModel(model,path):
    model.save(path)

def LoadModel(path):
    return keras.models.load_model(path)






