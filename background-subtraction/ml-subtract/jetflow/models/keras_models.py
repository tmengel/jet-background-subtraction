from __future__ import absolute_import

try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    
except ImportError as e:
    print('Could not import - ' + str(e))

__all__ = ['dnn', 'snn']

def dnn(input_features):
    inputlen = len(input_features)
    return keras.Sequential(
            [
                layers.Dense(100,input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(100, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(50, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(1, activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
            ])
    
def snn(input_feature):
    inputlen = len(input_feature)
    return keras.Sequential(
            [
                layers.Dense(inputlen+1,input_shape=(inputlen,),activation='relu',kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)),
                layers.Dense(1, activation='linear',kernel_regularizer=keras.regularizers.l2(0.001))
            ])

    
