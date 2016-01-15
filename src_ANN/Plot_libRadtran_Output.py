print '############################################################################################'
print '#                                                                                          #'
print '#  Plot_libRadtran_Output ( out_file, out_file_abs=None, fig_title="libRadtran output" )   #'
print '#                                                                                          #'
print '#  out_file     <str>        ::  filepath plus filename of libRadtran "flx"-output file    #'
print '#  out_file_abs <str> (opt)  ::  filepath plus filename of libRadtran "abs"-output file    #'
print '#  fig_title    <str> (opt)  ::  title of output figure                                    #'
print '#                                                                                          #'
print '#------------------------------------------------------------------------------------------#'
print '#                                                                                          #'
print '#  Compare_libRadtran_Output ( ANN_out, LUT_Out, fig_title=None )                          #'
print '#                                                                                          #'
print '#  ANN_out   <str>           ::  filepath plus filename of ANN "flx"-output file           #'
print '#  LUT_out   <str>           ::  filepath plus filename of LUT "flx"-output file           #'
print '#  fig_title <str> (opt)     ::  title of output figure                                    #'
print '#                                                                                          #'
print '############################################################################################'


def RMSE ( arr_1, arr_2, weights=None ):

    import numpy as np
  
    if np.shape(arr_1)!=np.shape(arr_2): raise ValueError('Arrays do not have same shapes :: {} .ne. {}'.format(np.shape(arr_1),np.shape(arr_2)) )
    s = np.sqrt ( np.average( (arr_1-arr_2)**2, weights=weights ) )
  
    return np.array([s, s/np.maximum(1.0E-8, abs(np.average(arr_2, weights=weights)))*100.0])



def Getting_Flux_Fields ( out_file ):

    import numpy as np
    
    data = np.loadtxt( out_file )
  
    _, x, y, z, edir, edn, eup, _, _, _ = data.T
    edir, edn, eup = [dfield.reshape( (np.max(z)+1,np.max(y)+1,np.max(x)+1) ) for dfield in edir, edn, eup]
    
    return edir, edn, eup



def Getting_Abso_Field ( out_file ):
    
    import numpy as np
  
    data = np.loadtxt( out_file )
  
    _, x, y, z, abso = data.T
    abso = np.reshape(abso, (np.max(z)+1,np.max(y)+1,np.max(x)+1) )
  
    return abso



def Plot_libRadtran_Output ( out_file, out_file_abs=None, fig_title=None ):
  
    import numpy as np
    import matplotlib.pyplot as plt
    
    edir, edn, eup = Getting_Flux_Fields ( out_file )
    rad = [edir,edn,eup]; var_name=['edir','edn','eup']
  
    if out_file_abs!=None:
        rad.append( Getting_Abso_Field(out_file_abs) )
        var_name.append('abso')
  
    fig = plt.figure()
    for num in range(len(rad)):
        ax  = fig.add_subplot(1, len(rad), num+1)
        img = ax.imshow(rad[num][:,:,0], origin='lower')
        plt.title( var_name[num] )
        plt.colorbar(img, fraction=0.046, pad=0.04)
    
    if fig_title!=None: fig.suptitle(fig_title, fontsize=18)
    plt.show()
  
    return



def Min_Max ( field_1, field_2 ):

    import numpy as np
  
    vmin = np.minimum( np.min(field_1), np.min(field_2) )
    vmax = np.maximum( np.max(field_1), np.max(field_2) )
  
    return vmin, vmax



def Error_Output ( rad ):

    print '######################'
    print '----------------------'
    var_name = ['edir','edn','eup','abso']
    for out in range(len(rad[0])):
        err, rel_err = RMSE (rad[0][out], rad[1][out])
        print '{} error'           .format(var_name[out])
        print 'abs err    ::    {}'.format(err)
        print 'rel err    ::    {}'.format(rel_err)
        print '----------------------'
    
    print '######################'
  
    return



def Compare_libRadtran_Output ( ANN_out, LUT_out, fig_title=None ):
  
    import numpy as np
    import matplotlib.pyplot as plt
  
    out_files = [ANN_out, LUT_out]
    rad=[]; var_name = ['edir','edn','eup','abs']; run = ['ANN','LUT']
    
    for out in range(len(out_files)):
        edir, edn, eup = Getting_Flux_Fields ( out_files[out] )
        rad.append( [edir, edn, eup] )
        try:
            out_file_abs = out_files[out].replace('flx','abs')
            abso = Getting_Abso_Field(out_file_abs)
            rad[out].append ( abso )
        except:
            print 'no absorption file available'
  
    vmin, vmax = [], []; ii=len(rad); jj=len(rad[0])
    for i in range( jj ):
        vminx, vmaxx = Min_Max ( rad[0][i], rad[1][i] )
        vmin.append(vminx); vmax.append(vmaxx)
    
    fig = plt.figure(); num=1
    for i in range( ii ):
        for j in range( jj ):
            ax  = fig.add_subplot(ii,jj,num); num=num+1
            img = ax.imshow(rad[i][j][:,:,0], origin='lower', vmin=vmin[j], vmax=vmax[j])
            plt.colorbar(img, fraction=0.046, pad=0.04)
            plt.title(run[i] + ' :: ' + var_name[j])
  
    if fig_title!=None: fig.suptitle(fig_title, fontsize=18)
    plt.show()
  
    Error_Output ( rad )
  
    return
