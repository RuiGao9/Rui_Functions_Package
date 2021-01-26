def MLP_Definer(input_dim=4,
                LayerNumber=1,
                neurons_1=4,activation_1='relu',
                neurons_2=4,activation_2='relu',
                neurons_3=4,activation_3='relu',
                neurons_4=4,activation_4='relu',
                neurons_out=2,activation_out='relu',
                lossfun='mse', opti='sgd'):
    '''
    The maximum hidden layers are 4 (4-layer total), and this MLP definer is designed for regression problem.
    parameter input_dim:
    parameter LyaerNumber: how many 
    
    '''
    import tensorflow as tf
    from keras.models import Sequential
    from keras.layers import Dense
    import numpy as np
    import numpy
    
    if LayerNumber==1:
        model = Sequential()
        model.add(Dense(neurons_1, input_dim=input_dim, activation=activation_1))
        model.add(Dense(neurons_out, activation=activation_out))
        model.compile(loss=lossfun, optimizer=opti, metrics=[tf.keras.metrics.MeanAbsoluteError()])
    elif LayerNumber==2:
        model = Sequential()
        model.add(Dense(neurons_1, input_dim=input_dim, activation=activation_1))
        model.add(Dense(neurons_2, activation=activation_2))
        model.add(Dense(neurons_out, activation=activation_out))
        model.compile(loss=lossfun, optimizer=opti, metrics=[tf.keras.metrics.MeanAbsoluteError()])
    elif LayerNumber==3:
        model = Sequential()
        model.add(Dense(neurons_1, input_dim=input_dim, activation=activation_1))
        model.add(Dense(neurons_2, activation=activation_2))
        model.add(Dense(neurons_3, activation=activation_3))
        model.add(Dense(neurons_out, activation=activation_out))
        model.compile(loss=lossfun, optimizer=opti, metrics=[tf.keras.metrics.MeanAbsoluteError()])
    elif LayerNumber==4:
        model = Sequential()
        model.add(Dense(neurons_1, input_dim=input_dim, activation=activation_1))
        model.add(Dense(neurons_2, activation=activation_2))
        model.add(Dense(neurons_3, activation=activation_3))
        model.add(Dense(neurons_4, activation=activation_4))
        model.add(Dense(neurons_out, activation=activation_out))
        model.compile(loss=lossfun, optimizer=opti, metrics=[tf.keras.metrics.MeanAbsoluteError()])
    else:
        print("Please make sure about the structure of the MLP model.")
    
    return(model)
    