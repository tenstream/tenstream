from functions import *
import numpy as np
import ffnet as ff
import matplotlib.pyplot as plt



FilePath    = '/usr/users/max.eckl/Dokumente/Arbeit/'
ArrayFolder = 'arrays/'
NNFolder    = 'networks/'

InputFile   = 'shuffled_input_array.npy'
TargetFile  = 'shuffled_target_array.npy'
NNName      = 'network.net'

train_percent = 0.1
train_steps   = 0.01
hidden_layers = 100







input_array  = np.load ( FilePath + ArrayFolder + InputFile  )
target_array = np.load ( FilePath + ArrayFolder + TargetFile )



net = Neural_Network(input_array, target_array, hidden_layers, train_percent, train_steps)


ff.savenet(net, FilePath + NNFolder + NNName)