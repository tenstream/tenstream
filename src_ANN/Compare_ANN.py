print '##################################################################'
print '#                                                                #'
print '#  Compare_ANN ( net_file, src, trg )                            #'
print '#                                                                #'
print '#  net_file <list>      ::  list of files which were saved with  #'
print '#                           "ffnet.savenet( net_file_2 )"        #'
print '#  src <numpy.ndarray>  ::  input array for ANN                  #'
print '#  trg <numpy.ndarray>  ::  target_array (LUT) for ANN           #'
print '#                                                                #'
print '##################################################################'


def Compare_ANN ( net_file, src, trg ):

  import numpy as np
  import ffnet as ff
  from helper_func import RMSE

  
  print '####################################'
  print '------------------------------------'

  net = []; net_out = []; rmse = []; rel_rmse = []
  for network in range(len(net_file)):
    net    .append( ff.loadnet( net_file[network] ) )

    out, _  = net[network].test(src, trg, iprint=0)
    net_out.append( out )
  
    rrmse, rrel_rmse = RMSE ( net[network](src), trg )
    rmse.append(rrmse); rel_rmse.append(rrel_rmse)


    print '                                     '
    print ' network_{}  ::  net_name = {}       '.format(network+1, net_file[network])
    print '                    RMSE = {}        '.format(rmse[network])
    print '                rel_RMSE = {}        '.format(rel_rmse[network])
    print '                                     '
    print '-------------------------------------'

  print '#####################################'

  return
