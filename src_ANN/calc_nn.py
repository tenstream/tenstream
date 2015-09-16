from ANN_functions import *
import numpy as np
import ffnet as ff




def Calc_Network_Diffuse ( hidden_neurons, train_percent, trainings ) :
    
    LUT_name = 'LUT_dstorder_8_10.diffuse.dx50.pspace.dz20.kabs20.ksca20.g3.delta_T_1.000.nc'
    Getting_Arrays_Diffuse ( LUT_name, dx=50, dy=50, LUT_path='LUT/', output_path='arrays/diffuse/')

    src = np.load('arrays/diffuse/s_src.npy')
    S   = np.load('arrays/diffuse/s_S.npy'  )
    
    net = Neural_Network(src, S, hidden_neurons)    
    
    Train_Network(src, S, net, trainings=trainings, train_percent=train_percent,
                  directory='networks/diffuse/', net_name='net_diffuse_', log_name='diffuse',
                  hidden_neurons=hidden_neurons, information='diffuse|dx 50|dy 50|\nLUT: %s' %('LUT/'+LUT_name))
    
    return 0


def Calc_Network_Direct ( hidden_neurons_S, hidden_neurons_T, train_percent, trainings ):
    
    LUT_name = 'LUT_dstorder_8_10.direct.dx50.pspace.dz20.kabs20.ksca20.g3.phi10.theta19.delta_T_1.000.nc'
    Getting_Arrays_Direct ( LUT_name, dx=50, dy=50, phi=0, theta=0, LUT_path='LUT/', output_path='arrays/direct/')
    
    src = np.load('arrays/direct/s_src.npy')
    S   = np.load('arrays/direct/s_S.npy'  )
    T   = np.load('arrays/direct/s_T.npy'  )
    
    net_S = Neural_Network(src, S, hidden_neurons_S)
    net_T = Neural_Network(src, T, hidden_neurons_T)
    
    Train_Network(src, S, net_S, trainings=trainings, train_percent=train_percent,
                  directory='networks/direct/', net_name='net_direct_S', log_name='direct_S',
                  hidden_neurons=hidden_neurons_S, information='S|direct|dx 50|dy 50|phi 0|theta 0\nLUT: %s' %('LUT/'+LUT_name))
    
    Train_Network(src, T, net_T, trainings=trainings, train_percent=train_percent,
                  directory='networks/direct/', net_name='net_direct_T', log_name='direct_T',
                  hidden_neurons=hidden_neurons_T, information='T|direct|dx 50|dy 50|phi 0|theta 0\nLUT: %s' %('LUT/'+LUT_name))
                  
    return 0

Calc_Network_Diffuse(hidden_neurons=65, train_percent=0.6, trainings=5)
Calc_Network_Direct(hidden_neurons_S=50, hidden_neurons_T=45, train_percent=0.6, trainings=5)
