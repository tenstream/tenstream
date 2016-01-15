print '####################################################################################'
print '#                                                                                  #'
print '#  Plot_Scatter_ANN_LUT ( net_file, src, trg, coeff )                              #'
print '#                                                                                  #'
print '#  net_file <str>       ::  file which was saved with "ffnet.savenet( net_file )"  #'
print '#  src <numpy.ndarray>  ::  input array for ANN                                    #'
print '#  trg <numpy.ndarray   ::  target array (LUT) for ANN                             #'
print '#  coeff <int>          ::  coefficient which will be plotted                      #'
print '#                                                                                  #'
print '####################################################################################'


def RMSE ( arr_1, arr_2, weights=None ):

    import numpy as np

    if np.shape(arr_1)!=np.shape(arr_2): raise ValueError('Arrays do not have same shapes :: {} .ne. {}'.format(np.shape(arr_1),np.shape(arr_2)))
    s = np.sqrt ( np.average( (arr_1-arr_2)**2, weights=weights ) )

    return np.array([s, s/np.maximum(1.0E-8, abs(np.average(arr_2, weights=weights)))*100.0])


def Getting_Arrays ( LUT_file ):
    
    import netCDF4 as nc
    import numpy as np

    LUT_data = nc.Dataset ( LUT_file, 'r' )
    
    coeff_type = 'diffuse'
    var        = ['dz', 'kabs', 'ksca', 'g']
    res        = Get_Dim( LUT_file, 'dx', 'dy', *var )
    LUT_names  = ['S', 'S_tol']
    LUT_prefix = 'dx{0:d}.dy{1:d}'.format(res['dx'],res['dy'])
    index, src, S, S_tol = [], [], [], []

    if LUT_file.rfind('direct')>=0:
        coeff_type = 'direct'
        res.update(Get_Dim( LUT_file, 'phi', 'theta' ))
        LUT_names.extend(['T','T_tol'])
        LUT_prefix = LUT_prefix+'.phi{0:d}.theta{1:d}'.format(res['phi'],res['theta'])
        T, T_tol = [], []
    
    data = [LUT_data.variables['{0:s}.dx{1:d}.dy{2:d}.pspace.{3:s}'.format(coeff_type,res['dx'],res['dy'],xvar)][:] for xvar in var]
    LUT  = [LUT_data.variables['{0:s}.{1:s}.{2:s}'.format(coeff_type,LUT_prefix,LUT_name)][:] for LUT_name in LUT_names]
    
    for  dz  in range( res[ 'dz' ] ):
        for kabs in range( res['kabs'] ):
            for ksca in range( res['ksca'] ):
                for  g   in range( res[ 'g'  ] ):

                    index.append ( [dz,kabs,ksca,g] )
                    src  .append ( [data[0][dz], data[1][kabs], data[2][ksca], data[3][g]] )
                    S    .append ( LUT[0][g,ksca,kabs,dz] )
                    S_tol.append ( LUT[1][g,ksca,kabs,dz] )
                    if coeff_type=='direct':
                        T    .append ( LUT[2][g,ksca,kabs,dz] )
                        T_tol.append ( LUT[3][g,ksca,kabs,dz] )

    dic = {'index':np.array(index), 'S':np.array(S), 'S_tol':np.array(S_tol)}; src=np.array(src)
    for ind, xvar in enumerate(var): dic[xvar]=np.array(data[ind])
    if coeff_type=='direct': 
        dic.update( {'T':np.array(T), 'T_tol':np.array(T_tol)} )
        app = np.append ( np.ones((len(src),1),dtype=float)*res['phi'], np.ones((len(src),1),dtype=float)*res['theta'], axis=1 )
        src = np.append ( src, app, axis=1 )
    dic.update( {'src':src} )

    return dic
  

def Plot_Scatter_ANN_LUT ( net_file, LUT_file, coeff, dir2dir=False ):

    import numpy as np
    import matplotlib.pyplot as plt
    import ffnet as ff

    LUT = Getting_Arrays(LUT_file)
    net = ff.loadnet ( net_file )

    src = LUT['src']
    trg = LUT['T'] if dir2dir else LUT['S']

    net_out = net(src)
    rel_error = np.minimum(np.abs((net_out[:,coeff]-trg[:,coeff])/trg[:,coeff]), 1.1)*100.0

    rmse, rel_rmse = RMSE(net_out[:,coeff], trg[:,coeff])

    plt.plot([-0.1,1.1],[-0.1,1.1], color='black')
    plt.scatter(net_out[:,coeff], trg[:,coeff],s=3, cmap='RdYlGn_r', c=rel_error, lw=0.0)
    ticks=np.linspace(0.0,100.0,11)
    cb = plt.colorbar(ticks=ticks); cb.set_label('rel. error [%]')
      
    plt.xlim([-0.1,1.1]); plt.ylim([-0.1,1.1]) 
    plt.title('%d    ::   %f   ::   %f' %(coeff,  rmse, rel_rmse))
    plt.grid()
    
    
    return
