import netCDF4 as nc
import numpy as np
import os as os
import ffnet as ff
import matplotlib.pyplot as plt



def Export_NN_to_NetCDF ( network, filename='NN.nc', filepath='/usr/users/max.eckl/Dokumente/Arbeit/networks/' ):
    
    FilePath = filepath
    FileName = filename
    File     = filepath + filename
    
    if not os.path.isdir(FilePath):
        os.mkdir(FilePath)
    if os.path.isfile(File):
        os.remove(File)
    
    dataset = nc.Dataset(File, 'w', format='NETCDF4')
    
    dataset.createDimension ( 'weights'  )
    dataset.createDimension ( 'conec'    ); dataset.createDimension ( 'conec2'    )
    dataset.createDimension ( 'units'    )
    dataset.createDimension ( 'inno'     )
    dataset.createDimension ( 'outno'    )
    dataset.createDimension ( 'eni'      ); dataset.createDimension ( 'eni2'      )
    dataset.createDimension ( 'deo'      ); dataset.createDimension ( 'deo2'      )
    dataset.createDimension ( 'inlimits' ); dataset.createDimension ( 'inlimits2' )
    
    weights  = dataset.createVariable('weights' , 'f8',  'weights'              )
    conec    = dataset.createVariable('conec'   , 'i' , ('conec', 'conec2'     ))
    units    = dataset.createVariable('units'   , 'f8',  'units'                )
    inno     = dataset.createVariable('inno'    , 'i' ,  'inno'                 )
    outno    = dataset.createVariable('outno'   , 'i' ,  'outno'                )
    eni      = dataset.createVariable('eni'     , 'f8', ('eni', 'eni2'         ))
    deo      = dataset.createVariable('deo'     , 'f8', ('deo', 'deo2'         ))
    inlimits = dataset.createVariable('inlimits', 'f8', ('inlimits','inlimits2'))
        
    weights [:] = network.weights
    conec   [:] = network.conec
    units   [:] = network.units
    inno    [:] = network.inno
    outno   [:] = network.outno
    eni     [:] = network.eni
    deo     [:] = network.deo
    inlimits[:] = network.inlimits

    dataset.close()
    
    
    return 0








def Getting_Input_Target ( LUT_filename, input_array_filename='input_array.npy',
                          target_array_filename='target_array.npy',
                          data_array_filename='data_array.npy',
                          output_filepath='arrays/',
                          LUT_filepath='/home/users/max.eckl/Dokumente/Arbeit/LUT/' ):


    FilePath = LUT_filepath
    FileName = LUT_filename
    File     = FilePath + FileName


    data       = nc.Dataset(File, 'r')
    data_array = data.variables['diffuse.dx40.dy40.S'][:]


    dim_g    = data_array.shape[0]
    dim_kabs = data_array.shape[1]
    dim_kref = data_array.shape[2]
    dim_dz   = data_array.shape[3]
    input_list  = []
    target_list = []

    for i in range(dim_g):
        for j in range(dim_kabs):
            for k in range(dim_kref):
                for l in range(dim_dz):
                    
                    input_list.append([i,j,k,l])
                    target_list.append(data_array[i,j,k,l].tolist())


    input_array  = np.array(input_list)
    target_array = np.array(target_list)
    
    
    if not os.path.isdir(output_filepath):
        os.mkdir(output_filepath)
    np.save(output_filepath+input_array_filename,  input_array )
    np.save(output_filepath+target_array_filename, target_array)
    np.save(output_filepath+data_array_filename,   data_array  )
    
    return 0
    



def Shuffle ( input_array, target_array ):
    
    indices = np.arange(0, input_array.shape[0])
    np.random.shuffle(indices)
    
    shuffled_input_array  = input_array[indices,:]
    shuffled_target_array = target_array[indices,:]
    
    np.save('arrays/shuffled_input_array', shuffled_input_array)
    np.save('arrays/shuffled_target_array', shuffled_target_array)
    
    return shuffled_input_array, shuffled_target_array




def RSE ( input_array_1, input_array_2 ):
    
    array_1 = np.array(input_array_1)
    array_2 = np.array(input_array_2)
    
    RSE = np.mean( np.abs( array_1 - array_2 ), axis=1 )
    
    return RSE






def RMSE ( input_array_1, input_array_2 ):
    
    array_1 = np.array(input_array_1)
    array_2 = np.array(input_array_2)
    
    RMSE = np.sqrt( np.mean ( (array_1 - array_2)**2 ) )
    
    return RMSE
    
    
    
    
    
    
    
    
    
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
       
    
    error = []
    error_2 = []
    err = 1.0E40
    j = 0
    
    
    while err >= RMSE(net(input_array[:(2*len_train)]), target_array[:(2*len_train)]):
        
        j = j+1
        fig = plt.figure(num=j)
        err = RMSE(net(input_array[:(2*len_train)]), target_array[:(2*len_train)])
        print 'RMSE ', err, ' for train session ', j
        
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