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