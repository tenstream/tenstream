from functions import *
import numpy as np
import ffnet as ff
import matplotlib.pyplot as plt



FilePath    = '/usr/users/max.eckl/Dokumente/Arbeit/'
ArrayFolder = 'arrays/'
NNFolder    = 'networks/'

InputFile   = 'shuffled_input_array.npy'
TargetFile  = 'shuffeld_target_array.npy'
NNName      = 'network.net'

train_percent = 0.1
train_steps   = 0.01
hidden_layers = 100



input_array  = np.load ( FilePath + ArrayFolder + InputFile  )
target_array = np.load ( FilePath + ArrayFolder + TargetFile )



net = Neural_Network(input_array, target_array, hidden_layers, train_percent, train_steps)


ff.savenet(net, FilePath + NNFolder + NNName)




def Neural_Network ( input_array, target_array , num_hidden_nodes, train_percent, train_steps, test_print=0 ):
    
    if type(input_array)==np.ndarray:
        input_array = input_array.tolist()
    
    if type(target_array)==np.ndarray:
        target_array = target_array.tolist()
    
    
    num_input_nodes    =    len (input_array [0])
    num_output_nodes   =    len (target_array[0])
    len_data           =    len (input_array    )
        
        
    len_train  = int(train_percent * len_data)
    num_steps  = int(train_percent / train_steps) 
    train_step = int(len_train/num_steps)
    
    conec = ff.mlgraph( (num_input_nodes, num_hidden_nodes, num_output_nodes) )
    net   = ff.ffnet( conec )
    
    fig = plt.figure()    
    
    error = []
    error_2 = []
    err = 1.0E40
    j = 0
    
    
    while err >= RMSE(net(input_array[:(2*len_train)]), target_array[:(2*len_train)]):
        
        j = j+1
        err = RMSE(net(input_array[:(2*len_train)]), target_array[:(2*len_train)])
        print 'RMSE ', err, ' for train session ' j
        
        for i in range(num_steps):
            net.train_tnc(input_array[(i*train_step):((i+1)*train_step)],target_array[(i*train_step):((i+1)*train_step)],nproc=8)
            print 'step ', i, ' is done!'
            
            error.append(RSE(net(input_array), target_array))
            error_2.append(RMSE(net(input_array), target_array))
            
            ax = fig.add_subplot(2,5,i+1)
            plt.plot(error[i])
            plt.title('step %d' %(i+1))
            
        plt.show()
        
    
    
    
    '''
    i=0
    error = 1.0
    while error >= 0.1:
        i=i+1
        net.train_tnc( input_array[:len_train], target_array[:len_train] )
        error = np.max(np.abs((net(input_array[len_train:(len_train+len_test)])-target_array[len_train:(len_train+len_test)])))
        error_index = np.argmax(np.abs((net(input_array[len_train:(len_train+len_test)])-target_array[len_train:(len_train+len_test)])))
        test = np.max(error/np.maximum(target_array[error_index],1.0E-30))
        print i, error, target_array[error_index], net(input_array[error_index]), test
        '''
    '''
    i=0
    error = -1.0  
    while RMSE(net(input_array), target_array) != error:
        i = i+1
        error = RMSE(net(input_array), target_array) 
        net.train_tnc( input_array, target_array, nproc=8 , maxfun=10000)
        if i >= 20:
            break
        
    
    output_array, regression = net.test(input_array, target_array, iprint=test_print)
    '''
    
    return net
