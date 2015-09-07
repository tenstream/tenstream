###############################################################################
#-----------------------------------------------------------------------------#
#---------functions for calculating, training, and testing of-----------------#
#-------------------artificial neural networks--------------------------------#
#-----------------------------------------------------------------------------#
###############################################################################


def Shuffle_2D_X ( array ):
    
    import numpy as np
    
    if type(array)==tuple or type(array)==list:
        indices = np.arange(0, array[0].shape[0],dtype=int)
        np.random.shuffle(indices)        
        
        s_array = []
        for a in range(len(array)):
            if len(array[a].shape)==2:
                s_array.append(array[a][indices,:])
            else:
                print '######################################'
                print '!!! ERROR IN FUNCTION Shuffle_2D_X !!!'
                print '---   array %d has wrong dimenson   ---' %a
                print '######################################'
            
        return s_array
        
    elif type(array)==np.ndarray:
        indices = np.arange(0, array.shape[0],dtype=int)
        np.random.shuffle(indices)
        
        if len(array.shape)==2:
            s_array = array[indices,:]
            return s_array
        else:
            print '######################################'
            print '!!! ERROR IN FUNCTION Shuffle_2D_X !!!'
            print '---    array has wrong dimenson    ---' 
            print '######################################'
            return 0
    
    else:
        print '#########################################################################'
        print '!!!                  ERROR IN FUNCTION Shuffle_2D_X                   !!!'
        print '--- input data type is %s but tuple, list or np.ndarray is needed ---' %str(type(array))
        print '#########################################################################'
        return 0
    



def RMSE ( array_1, array_2 ):
    
    import numpy as np
    
    
    if type(array_1)==list:
        array_1 = np.array(array_1)
    if type(array_2)==list:
        array_2 = np.array(array_2)
        
    if array_1.shape!=array_2.shape:
        print '#########################################'
        print '!!!       ERROR IN FUNCTION RMSE      !!!'
        print '--- shapes of arrays are not the same ---'
        print '#########################################'
        return 0
    
    RMSE = np.sqrt( np.mean( (array_1 - array_2)**2 ) )
    
    
    return RMSE

def rmse(a,b,weights=None):
        import numpy as np
        if np.shape(a)!=np.shape(b): raise ValueError('Arrays do not have same shapes :: {} .ne. {}'.format(np.shape(a),np.shape(b)) )

        s= np.sqrt ( np.average( (a - b)**2, weights=weights ) )
        return np.array([s,s/np.maximum(1e-8,abs(np.average(b,weights=weights)))*100])





def Export_NN_to_NetCDF ( network, FileName='neural_network.nc', FilePath='' ):
    
    import netCDF4 as nc
    import numpy as np
    import os as os
    
    
    File = FilePath +'/'+ FileName
    print 'Exporting Network to .nc file :: ',File
    
    if not os.path.isdir(FilePath):
        os.mkdir(FilePath)
    if os.path.isfile(File):
        os.remove(File)
    
    dataset = nc.Dataset(File, 'w', format='NETCDF4')
    
    dataset.createDimension ( 'weights.dim1'  )
    dataset.createDimension ( 'conec.dim1'    ); dataset.createDimension ( 'conec.dim2'    )
    dataset.createDimension ( 'units.dim1'    )
    dataset.createDimension ( 'inno.dim1'     )
    dataset.createDimension ( 'outno.dim1'    )
    dataset.createDimension ( 'eni.dim1'      ); dataset.createDimension ( 'eni.dim2'      )
    dataset.createDimension ( 'deo.dim1'      ); dataset.createDimension ( 'deo.dim2'      )
    dataset.createDimension ( 'inlimits.dim1' ); dataset.createDimension ( 'inlimits.dim2' )
    
    weights  = dataset.createVariable('weights' , 'f8',  'weights.dim1'                  )
    conec    = dataset.createVariable('conec'   , 'i' , ('conec.dim1'   , 'conec.dim2'  ))
    units    = dataset.createVariable('units'   , 'f8',  'units.dim1'                    )
    inno     = dataset.createVariable('inno'    , 'i' ,  'inno.dim1'                     )
    outno    = dataset.createVariable('outno'   , 'i' ,  'outno.dim1'                    )
    eni      = dataset.createVariable('eni'     , 'f8', ('eni.dim1'     , 'eni.dim2'    ))
    deo      = dataset.createVariable('deo'     , 'f8', ('deo.dim1'     , 'deo.dim2'    ))
    inlimits = dataset.createVariable('inlimits', 'f8', ('inlimits.dim1','inlimits.dim2'))
        
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




def Getting_Arrays_Diffuse ( LUT_name, dx, dy, LUT_path='', index_name='index', 
                             src_name='src', S_name='S', output_path='arrays/',
                             messages=True ):
    
    import netCDF4 as nc
    import numpy as np
    import os as os


    File     = LUT_path + LUT_name

    LUT_data = nc.Dataset(File, 'r')
    LUT      = LUT_data.variables['diffuse.dx%d.dy%d.S' %(dx,dy)][:]

    var_name = ['g', 'ksca', 'kabs', 'dz']; data = []
    for var in var_name:
        data.append( LUT_data.variables ['diffuse.dx%d.dy%d.pspace.%s'  %(dx,dy,var)][:] )
    
    
    index_list = []; src_list = []; S_list = []
    for g in range(LUT.shape[0]):
        for ksca in range(LUT.shape[1]):
            for kabs in range(LUT.shape[2]):
                for dz in range(LUT.shape[3]):
                    
                    index_list.append([g,ksca,kabs,dz])
                    src_list  .append([data[0][g], data[1][ksca], data[2][kabs], data[3][dz]])
                    S_list    .append(LUT[g,ksca,kabs,dz])

    index = np.array( index_list )
    src   = np.array( src_list   )
    S     = np.array( S_list     )
    
    
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    np.save( output_path + index_name  + '.npy', index )
    np.save( output_path + src_name    + '.npy', src   )
    np.save( output_path + S_name      + '.npy', S     )
    
    
    s_arrays = Shuffle_2D_X( [index, src, S] )
    np.save( output_path + 's_' + index_name  + '.npy', s_arrays[0] )
    np.save( output_path + 's_' + src_name    + '.npy', s_arrays[1] )
    np.save( output_path + 's_' + S_name      + '.npy', s_arrays[2] )    
    
    
    if messages:
        print '######################################################'
        print 'Getting_Arrays_Diffuse:'
        print 'LUT filename\t\t%s'                               %LUT_name
        if LUT_path=='':
            print 'LUT filepath\t\t%s'                           %('THIS_DIRECTORY')
        else:
            print 'LUT filepath\t\t%s'                           %LUT_path
        print '(dx, dy)\t\t\t(%d, %d)'                           %(dx, dy)
        print 'output filepath\t\t%s'                            %output_path
        print 'output arrays\t\t%s.npy, %s.npy, %s.npy'          %(index_name, src_name, S_name)
        print 'shuffled output arrays\t\t%s.npy, %s.npy, %s.npy' %('s_'+index_name, 's_'+src_name, 's_'+S_name)
        print 'dimension of index arrays\t(%d, %d)'              %(index.shape[0], index.shape[1])
        print 'dimension of src arrays\t(%d, %d)'                %(src.shape[0], src.shape[1])
        print 'dimension of S arrays\t\t(%d, %d)'                %(S.shape[0], S.shape[1])
        print 'variables\t\t\t%s, %s, %s, %s'                    %(var_name[0], var_name[1], var_name[2], var_name[3])
        print '######################################################'
        
    
    return 0




def Getting_Arrays_Direct ( LUT_name, dx, dy, phi, theta, LUT_path='', index_name='index',
                            src_name='src', S_name='S', T_name='T', output_path='arrays/',
                            messages=True ):
                           
    import netCDF4 as nc
    import numpy as np
    import os as os
    
    
    File = LUT_path + LUT_name
    
    LUT_data = nc.Dataset(File, 'r')
    LUT_T    = LUT_data.variables['direct.dx%d.dy%d.phi%d.theta%d.T' %(dx,dy,phi,theta)][:]
    LUT_S    = LUT_data.variables['direct.dx%d.dy%d.phi%d.theta%d.S' %(dx,dy,phi,theta)][:]
    
    
    var_name = ['g', 'ksca', 'kabs', 'dz']; data = []
    for var in var_name:
        data.append( LUT_data.variables ['direct.dx%d.dy%d.pspace.%s' %(dx,dy,var)][:] )
        
    
    index_list = []; src_list = []; T_list = []; S_list = []
    for g in range(LUT_T.shape[0]):
        for ksca in range(LUT_T.shape[1]):
            for kabs in range(LUT_T.shape[2]):
                for dz in range(LUT_T.shape[3]):
                    
                    index_list.append([g,ksca,kabs,dz])
                    src_list  .append([data[0][g], data[1][ksca], data[2][kabs], data[3][dz]])
                    T_list    .append(LUT_T[g,ksca,kabs,dz])
                    S_list    .append(LUT_S[g,ksca,kabs,dz])
    
    index = np.array( index_list )
    src   = np.array( src_list   )
    T     = np.array( T_list     )
    S     = np.array( S_list     )
    
    
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    np.save( output_path + index_name + '.npy', index )
    np.save( output_path + src_name   + '.npy', src   )
    np.save( output_path + T_name     + '.npy', T     )
    np.save( output_path + S_name     + '.npy', S     )
    
    
    s_arrays = Shuffle_2D_X( [index, src, T, S] )
    np.save( output_path + 's_' + index_name + '.npy', s_arrays[0] )
    np.save( output_path + 's_' + src_name   + '.npy', s_arrays[1] )
    np.save( output_path + 's_' + T_name     + '.npy', s_arrays[2] )
    np.save( output_path + 's_' + S_name     + '.npy', s_arrays[3] )
    
    
    if messages:
        print '######################################################'
        print 'Getting_Arrays_Direct:'
        print 'LUT filename\t\t%s'                                       %LUT_name
        if LUT_path=='':
            print 'LUT filepath\t\t%s'                                   %('THIS_DIRECTORY')
        else:
            print 'LUT filepath\t\t%s'                                   %LUT_path
        print '(dx, dy, phi, theta)\t\t(%d, %d, %d, %d)'                 %(dx, dy, phi, theta)
        print 'output filepath\t\t%s'                                    %output_path
        print 'output arrays\t\t%s.npy, %s.npy, %s.npy, %s.npy'          %(index_name, src_name, T_name, S_name)
        print 'shuffled output arrays\t\t%s.npy, %s.npy, %s.npy, %s.npy' %('s_'+index_name, 's_'+src_name, 's_'+T_name, 's_'+S_name)
        print 'dimension of index arrays\t(%d, %d)'                      %(index.shape[0], index.shape[1])
        print 'dimension of src arrays\t(%d, %d)'                        %(src.shape[0], src.shape[1])
        print 'dimension of T arrays\t\t(%d, %d)'                        %(T.shape[0], T.shape[1])
        print 'dimension of S arrays\t\t(%d, %d)'                        %(S.shape[0], S.shape[1])
        print 'variables\t\t\t%s, %s, %s, %s'                            %(var_name[0], var_name[1], var_name[2], var_name[3])
        print '######################################################'
            
    
    return 0




def Neural_Network ( src, trg, hidden_neurons, messages=True ):
    
    import ffnet as ff
    import numpy as np
    
    
    if messages:
        print '###########################################################'
        print '              calculating neural network...                '
    
    if type(src)==list:
        src_neurons    =    len ( src [0] )
    else:
        src_neurons    =    src.shape[1]    
    if type(src)==list:
        trg_neurons   =    len ( trg [0] )
    else:
        trg_neurons   =    trg.shape[1]
        
    if hidden_neurons<=min(src_neurons, trg_neurons):
        print '!!!                      WARNING                        !!!'
        print '--- too little hidden neurons; underfitting is probable ---' 
    if hidden_neurons>=max(src_neurons, trg_neurons):
        print '!!!                      WARNING                        !!!'
        print '---  too much hidden neurons; overfitting is probable   ---'          
              
    conec = ff.mlgraph( (src_neurons, hidden_neurons, trg_neurons) )
    net   = ff.ffnet( conec )
    
    if messages:
        print '           neural network has shape (%d,%d,%d)             '  %(src_neurons, hidden_neurons, trg_neurons)
        print '###########################################################'
    
    return net



 

def Train_Network ( src, trg, net, trainings, train_percent, nproc=4,
                    messages=True, save_log=True, save_network=True, directory='',
                    net_name='NOT_SET', log_name='NOT_SET', train_type='partly',
                    train_start=0, information='', hidden_neurons='UNKNOWN',
                    shuffle='TRUE' ):
                        
    import numpy as np
    import ffnet as ff
    import os as os
    

    if train_type!='complete' and train_type!='partly':
        print '################################################################'
        print '!!!                    ERROR IN Train_Network                !!!'
        print '--- unknown training type "%s"; use "complete" or "partly" ---' %train_type
        print '################################################################'
        return 0
    

    if type(src)==np.ndarray:
        src = src.tolist()
    if type(trg)==np.ndarray:
        trg = trg.tolist()
    
    if directory=='':
        data_dir = 'THIS_DIRECTORY'
    else:
        data_dir = directory
        if not os.path.isdir(directory):
            os.mkdir(directory)
        
    
    if log_name=='NOT_SET':
        log_name = 'log'
    if net_name=='NOT_SET':
        net_name = 'network'
    
    
    len_data    = len ( src )
    train_step  = int(train_percent * len_data)
    train_parts = int(1.0/train_percent)

    
    if messages:
        print '-------------------------------------------'
        print '        TRAINING NEURAL NETWORK            '
        print '-------------------------------------------'
        print '###########################################'
        print 'neural network configuration:'
        print 'hidden neurons\t\t%s'               %str(hidden_neurons)      
        print 'length of data\t\t%d'               %len_data
        print 'train step percent\t\t%f'           %train_percent
        print 'train step\t\t\t%d'                 %train_step
        if train_type=='complete':
            print 'train parts\t\t\t%d'            %train_parts
        print 'train type\t\t\t%s'                 %train_type
        if train_type=='complete':
            print 'trainings\t\t\t%d'              %trainings
        if train_type=='partly':
            print 'train start (index)\t\t%d'      %train_start
        print 'shape input array\t\t(%d, %d)'      %(len(src),len(src[0]))
        print 'shape target array\t\t(%d, %d)'     %(len(trg),len(trg[0]))
        print 'used shuffled arrays\t\t%s'         %shuffle
        print 'used processors\t\t%d'              %nproc
        print '- - - - - - - - - - - - - - - - - - - - - -'
        print 'log data name:\t\t%s'               %log_name+'.dat'
        print 'save directory:\t\t%s'              %data_dir
        print '###########################################\n'
        if information!='':
            print 'INFORMATION: %s\n'              %information
            print '###########################################\n'
        
        
        
    if save_log:
        log = open(directory + log_name + '.dat', 'w')
        log.write('###########################################\n'                   )
        log.write('neural network configuration:\n'                                 )
        log.write('hidden neurons\t\t\t%s\n'               %str(hidden_neurons)     )
        log.write('length of data\t\t\t%d\n'               %len_data                )
        log.write('train step percent\t\t%f\n'             %train_percent           )
        log.write('train step\t\t\t%d\n'                   %train_step              )
        if train_type=='complete':
            log.write('train parts\t\t\t%d\n'              %train_parts             )
        log.write('train type\t\t\t%s\n'                   %train_type              )
        if train_type=='complete':
            log.write('trainings\t\t\t%d\n'                %trainings               )
        if train_type=='partly':
            log.write('train start (index)\t\t%d\n'        %train_start             )
        log.write('shape input array\t\t(%d, %d)\n'        %(len(src),len(src[0]))  )
        log.write('shape target array\t\t(%d, %d)\n'       %(len(trg),len(trg[0]))  )
        log.write('used shuffled arrays\t\t%s\n'           %shuffle                 )
        log.write('used processors\t\t\t%d\n'              %nproc                   )
        log.write('- - - - - - - - - - - - - - - - - - - - - -\n'                   )
        log.write('log data name:\t\t\t%s\n'               %(log_name+'.dat')       )
        log.write('save directory:\t\t\t%s\n'              %data_dir                )
        log.write('###########################################\n\n'                 )
        if information!='':
            log.write('INFORMATION: %s\n\n'                %information             )
            log.write('###########################################\n\n'             )
    
    


    if train_type=='complete':
        
        for training in range(trainings):
            
            restart = 0
            
            if messages:
                print 'starting training session %d of %d' %(training+1, trainings)
                if save_network:
                    print 'training step\tfrom\tto\tRMSE_complete\tRMSE_trained\tRMSE_untrained\tsave_name'
                else:
                    print 'training step\tfrom\tto\tRMSE_complete\tRMSE_trained\tRMSE_untrained'
            if save_log:
                log.write('training session %d of %d\n' %(training+1, trainings))
                if save_network:
                    log.write('training step\tfrom\tto\tRMSE_complete\tRMSE_trained\tRMSE_untrained\tsave name\n')
                else:
                    log.write('training step\tfrom\tto\tRMSE_complete\tRMSE_trained\tRMSE_untrained\n')
            
            for part in range(train_parts):
                
                if train_start==0:
                    start =  part    * train_step
                    stop  = (part+1) * train_step
                else:
                    start = train_start +  part    * train_step - restart
                    stop  = train_start + (part+1) * train_step - restart
                    if start>=len(src):
                        restart = train_start + part*train_step
                        print restart
                        start   = 0
                        stop    = train_step
                if stop>len(src):
                    stop = len(src)

        
                net.train_tnc(src[start:stop], trg[start:stop], nproc=nproc)
            
                
                output_complete, regression_complete = net.test(src            , trg            , iprint=0)
                output_trained , regression_train    = net.test(src[start:stop], trg[start:stop], iprint=0)
                if start==0:
                    trg_untrained = trg[stop:]
                    output_untrained, regression_untrained = net.test(src[stop:], trg_untrained, iprint=0)
                elif stop==len(src):
                    trg_untrained = trg[:start]
                    output_untrained, regression_untrained = net.test(src[:start], trg_untrained, iprint=0)
                else:
                    src_untrained = np.append(np.array(src[:start]), np.array(src[stop:]), axis=0)
                    trg_untrained = np.append(np.array(trg[:start]), np.array(trg[stop:]), axis=0)
                    output_untrained, regression_untrained = net.test(src_untrained, trg_untrained, iprint=0)
                
                error_complete  = RMSE ( trg            , output_complete  )
                error_trained   = RMSE ( trg[start:stop], output_trained   )
                error_untrained = RMSE ( trg_untrained  , output_untrained )
        
                if save_network:
                    nsave = net_name + '%d.net' %(training+1)
                    ff.savenet(net, directory + nsave)
                else:
                    nsave = ''
                if messages:
                    print '%d\t\t%d\t%d\t%f\t\t%f\t\t%f\t\t%s' %((training+1), start, (stop-1), error_complete, error_trained, error_untrained, nsave)
                if save_log:
                    log.write('%d\t\t%d\t%d\t%f\t%f\t%f\t%s\n' %((training+1), start, (stop-1), error_complete, error_trained, error_untrained, nsave))
            
            if save_log:
                log.write('\n')
    
    
    
    if train_type=='partly':
            
        start = train_start
        if start+train_step<len(src):
            stop = start+train_step
        else:
            stop = len(src)
            print 'train range is too big!!! -> train_percent is set to %f instead of %f' %((stop-start)/len_data,train_percent)
            if save_log:
                log.write('!!!train range is too big!!!\n')
                log.write('train_percent is set to %f instead of %f\n' %((float(stop)-float(start))/float(len_data),train_percent))
                log.write('###########################################\n\n')
         
        if messages:
            print '%d trainings with input array from index %d to %d\n' %(trainings, start, stop)
            if save_network:
                print 'training\tRMSE_complete\tRMSE_trained\tRMSE_untrained\tsave name'
            else:
                print 'training\tRMSE_complete\tRMSE_trained\tRMSE_untrained'
        if save_log:
            log.write('%d trainings with input array from index %d to %d\n\n' %(trainings, start, stop))
            if save_network:
                log.write('training\tRMSE_complete\tRMSE_trained\tRMSE_untrained\tsave name\n')
            else:
                log.write('training\tRMSE_complete\tRMSE_trained\tRMSE_untrained\n')
                
        for training in range(trainings):
        
            net.train_tnc(src[start:stop], trg[start:stop], nproc=nproc)
            
            output_complete, regression_complete = net.test(src            , trg            , iprint=0)
            output_trained , regression_train    = net.test(src[start:stop], trg[start:stop], iprint=0)
            if start==0:
                trg_untrained = trg[stop:]
                output_untrained, regression_untrained = net.test(src[stop:], trg_untrained, iprint=0)
            elif stop==len(src):
                trg_untrained = trg[:start]
                output_untrained, regression_untrained = net.test(src[:start], trg_untrained, iprint=0)
            else:
                src_untrained = np.append(np.array(src[:start]), np.array(src[stop:]), axis=0)
                trg_untrained = np.append(np.array(trg[:start]), np.array(trg[stop:]), axis=0)
                output_untrained, regression_untrained = net.test(src_untrained, trg_untrained, iprint=0)
                
            error_complete  = RMSE ( trg            , output_complete  )
            error_trained   = RMSE ( trg[start:stop], output_trained   )
            error_untrained = RMSE ( trg_untrained  , output_untrained ) 
            
            if save_network:
                nsave = net_name + '%d.net' %(training+1)
                ff.savenet(net, directory + nsave)
            else:
                nsave = ''
            if messages:
                print '%d\t%f\t\t%f\t\t%f\t\t%s' %((training+1), error_complete, error_trained, error_untrained, nsave)
            if save_log:
                log.write('%d\t\t%f\t%f\t%f\t%s\n' %((training+1), error_complete, error_trained, error_untrained, nsave))
                
    
    if messages:
        print ''
        print '-------------------------------------------'
        print '  NETWORK HAS BEEN SUCCESSFULLY TRAINED    '
        print '-------------------------------------------'
    if save_log:
        log.close()
    
    
    return 0
