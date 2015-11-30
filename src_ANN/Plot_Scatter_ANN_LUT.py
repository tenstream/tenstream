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

  

def Plot_Scatter_ANN_LUT ( net, src, trg, coeff ):

    import numpy as np
    import matplotlib.pyplot as plt
    import ffnet as ff
    from helper_func import RMSE


    net = ff.loadnet ( net_file )

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
