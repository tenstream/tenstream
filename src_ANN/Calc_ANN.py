print '##############################################################################################################'
print '#                                                                                                            #'
print '#  Calc_ANN ( LUT_file, coeff_type, ANN_setup, ANN_basename="network", connectivity="standard",              #'
print '#             test_perc=0.3, err_inc=0.001, nproc=8, shuffle_alw=True, netcdf=None, info=None   )            #'
print '#                                                                                                            #'
print '##############################################################################################################'


def RMSE(a,b,weights=None):
  import numpy as np
  
  if np.shape(a)!=np.shape(b): raise ValueError('Arrays do not have same shapes :: {} .ne. {}'.format(np.shape(a),np.shape(b)) )
  s= np.sqrt ( np.average( (a - b)**2, weights=weights ) )
    
  return np.array([s,s/np.maximum(1e-8,abs(np.average(b,weights=weights)))*100])


def Shuffle_2D_X ( array ):
  import numpy as np
  
  if len(array.shape)!=2: raise ValueError('Array is not 2D :: dim(array) = {}'.format(len(array.shape)) )
  indices = np.arange(0, array.shape[0], dtype=int)
  np.random.shuffle(indices)
  s_array = array[indices,:]

  return indices, s_array


def Get_Dim ( LUT_file, key ):
  
  ind = LUT_file.rfind( key ) + len(key); num=0
  while LUT_file[ind+num]!='.': num=num+1
  dim = int(LUT_file[ind:ind+num])

  return dim


def Getting_Arrays_Diffuse ( LUT_file ):
  import netCDF4 as nc
  import numpy as np

  dx        = Get_Dim ( LUT_file, 'dx' )
  LUT_data  = nc.Dataset ( LUT_file, 'r' )
  LUT_S     = LUT_data.variables['diffuse.dx{}.dy{}.S'.format(dx, dx)][:]
  LUT_S_tol = LUT_data.variables['diffuse.dx{}.dy{}.S_tol'.format(dx, dx)][:]

  var_name=['dz', 'kabs', 'ksca', 'g']; data=[]
  for var in var_name:
    data.append( LUT_data.variables['diffuse.dx{}.dy{}.pspace.{}'.format(dx, dx, var)][:] )
  
  index, src, S, S_tol = [], [], [], []
  for dz   in range(LUT_S.shape[3]):
    for kabs in range(LUT_S.shape[2]):
      for ksca in range(LUT_S.shape[1]):
        for g   in range(LUT_S.shape[0]):
          index.append ( [dz,kabs,ksca,g] )
          src  .append ( [data[0][dz],data[1][kabs],data[2][ksca],data[3][g]] )
          S    .append ( LUT_S     [g,ksca,kabs,dz] )
          S_tol.append ( LUT_S_tol [g,ksca,kabs,dz] )
  
  index,src,S,S_tol = np.array(index),np.array(src),np.array(S),np.array(S_tol)

  return index, src, S, S_tol


def Getting_Arrays_Direct ( LUT_file ):
  import netCDF4 as nc
  import numpy as np

  dim=[]
  for dim_value in ['dx', 'dy', 'phi', 'theta']:
    dim.append( Get_Dim( LUT_file, dim_value ) )
  
  LUT_data = nc.Dataset(LUT_file, 'r')

  data=[]; LUT=[]
  data_name=['dz','kabs','ksca','g']; LUT_name=['S','S_tol','T','T_tol']
  for var in range(len(data_name)):
    data.append( LUT_data.variables['direct.dx{}.dy{}.pspace.{}'.format(dim[0],dim[1],data_name[var])][:] )
    LUT .append( LUT_data.variables['direct.dx{}.dy{}.phi{}.theta{}.{}'.format(dim[0],dim[1],dim[2],dim[3],LUT_name[var])][:] )

  index,src,S,S_tol,T,T_tol=[],[],[],[],[],[]
  for dz   in range( data[0].shape[0] ):
    for kabs in range( data[1].shape[0] ):
      for ksca in range( data[2].shape[0] ):
        for g    in range( data[3].shape[0] ):
          index.append ( [dz,kabs,ksca,g] )
          src  .append ( [data[0][dz],data[1][kabs],data[2][ksca],data[3][g],float(dim[2]),float(dim[3])] )
          S    .append ( LUT[0] [g,ksca,kabs,dz] )
          S_tol.append ( LUT[1] [g,ksca,kabs,dz] )
          T    .append ( LUT[2] [g,ksca,kabs,dz] )
          T_tol.append ( LUT[3] [g,ksca,kabs,dz] )

  index,src,S,S_tol,T,T_tol = np.array(index),np.array(src),np.array(S),np.array(S_tol),np.array(T),np.array(T_tol)

  return index, src, S, S_tol, T, T_tol


def Getting_Arrays ( LUT_file, coeff_type ):
  from numpy import append

  if type(LUT_file)==str:
    if coeff_type=='diff2diff' or coeff_type=='diffuse':
      _, src, trg, _ = Getting_Arrays_Diffuse ( LUT_file )
    elif coeff_type=='dir2diff':
      _, src, trg, _, _, _ = Getting_Arrays_Direct ( LUT_file )
    elif coeff_type=='dir2dir':
      _, src, _, _, trg, _ = Getting_Arrays_Direct ( LUT_file )
    else:
      raise ValueError('unknown coeff_type :: "{}"; use "diffuse", "diff2diff", "dir2diff", or "dir2dir"'.format(coeff_type))

  if type(LUT_file)==list:
    if coeff_type=='diff2diff' or coeff_type=='diffuse':
      _, src, trg, _ = Getting_Arrays_Diffuse ( LUT_file[0] )
      for LUT in range(len(LUT_file)-1):
        _, xsrc, xtrg, _ = Getting_Arrays_Diffuse ( LUT_file[LUT+1] )
        src, trg = append(src, xsrc, axis=0), append(trg, xtrg, axis=0)
    elif coeff_type=='dir2diff':
      _, src, trg, _, _, _ = Getting_Arrays_Direct ( LUT_file[0] )
      for LUT in range(len(LUT_file)-1):
        _, xsrc, xtrg, _, _, _ = Getting_Arrays_Direct ( LUT_file[LUT+1] )
        src, trg = append(src, xsrc, axis=0), append(trg, xtrg, axis=0)
    elif coeff_type=='dir2dir':
      _, src, _, _, trg, _ = Getting_Arrays_Direct ( LUT_file[0] )
      for LUT in range(len(LUT_file)-1):
        _, xsrc, _, _, xtrg, _ = Getting_Arrays_Direct ( LUT_file[LUT+1] )
        src, trg = append(src, xsrc, axis=0), append(trg, xtrg, axis=0)
    else:
      raise ValueError('unknown coeff_type :: "{}"; use "diffuse", "diff2diff", "dir2diff", or "dir2dir"'.format(coeff_type))
  
  return src, trg


def Init_ANN ( setup, connectivity='standard' ):
  import ffnet as ff

  if connectivity=='standard':
    conec = ff.mlgraph ( setup )
  elif connectivity=='full':
    conec = ff.tmlgraph( setup )
  else:
    raise ValueError('Unknown connectivity type "{}" :: available types are "standard" and "full"'.format(connectivity))
  net = ff.ffnet ( conec )

  return net


def Train_ANN ( net, src, trg, test_perc=0.3, nproc=8 ):

  ind = int ( src.shape[0]*(1.0-test_perc) )
  src_train, trg_train = src[:ind,:], trg[:ind,:]
  src_test , trg_test  = src[ind:,:], trg[ind:,:]

  src_list, trg_list = [src,src_train,src_test], [trg,trg_train,trg_test]
  
  if nproc=='ncpu':
    try:
      net.train_tnc ( src_train, trg_train, nproc='ncpu' )
    except:
      print 'nproc="ncpu" does not work; nproc=8 is used instead'
      net.train_tnc ( src_train, trg_train, nproc=8      )
  else:
    net.train_tnc ( src_train, trg_train, nproc=nproc  )

  err, rel_err = [], []
  for i in range(3):
    out, _ = net.test ( src_list[i], trg_list[i], iprint=0 )
    xerr, xrel_err = RMSE ( out, trg_list[i] )
    err.append(xerr); rel_err.append(xrel_err)
  
  return net, err, rel_err


def ANN_to_NetCDF ( net, out_file, iprint=True ):
  import netCDF4 as nc
  import numpy as np
  import ffnet as ff    
  
  if iprint: print 'Exporting Network "{}" to .nc file :: {}'.format(net, out_file)
  if type(net)==str:
    network = ff.loadnet ( net_file )
  else:
    network = net
  
  # transpose arrays to fortran order
  Tweights  = network.weights .T
  Tconec    = network.conec   .T
  Tunits    = network.units   .T
  Tinno     = network.inno    .T
  Toutno    = network.outno   .T
  Teni      = network.eni     .T
  Tdeo      = network.deo     .T
  Tinlimits = network.inlimits.T
  
  dataset = nc.Dataset(out_file, 'w', format='NETCDF4')
  
  dataset.createDimension ( 'weights_dim1' , np.shape(Tweights )[0] )
  dataset.createDimension ( 'conec_dim1'   , np.shape(Tconec   )[0] ); dataset.createDimension ( 'conec_dim2'   , np.shape(Tconec   )[1] )
  dataset.createDimension ( 'units_dim1'   , np.shape(Tunits   )[0] )
  dataset.createDimension ( 'inno_dim1'    , np.shape(Tinno    )[0] )
  dataset.createDimension ( 'outno_dim1'   , np.shape(Toutno   )[0] )
  dataset.createDimension ( 'eni_dim1'     , np.shape(Teni     )[0] ); dataset.createDimension ( 'eni_dim2'     , np.shape(Teni     )[1] )
  dataset.createDimension ( 'deo_dim1'     , np.shape(Tdeo     )[0] ); dataset.createDimension ( 'deo_dim2'     , np.shape(Tdeo     )[1] )
  dataset.createDimension ( 'inlimits_dim1', np.shape(Tinlimits)[0] ); dataset.createDimension ( 'inlimits_dim2', np.shape(Tinlimits)[1] )

  weights  = dataset.createVariable('weights' , 'f8',  'weights_dim1'                  )
  conec    = dataset.createVariable('conec'   , 'i' , ('conec_dim1'   , 'conec_dim2'  ))
  units    = dataset.createVariable('units'   , 'f8',  'units_dim1'                    )
  inno     = dataset.createVariable('inno'    , 'i' ,  'inno_dim1'                     )
  outno    = dataset.createVariable('outno'   , 'i' ,  'outno_dim1'                    )
  eni      = dataset.createVariable('eni'     , 'f8', ('eni_dim1'     , 'eni_dim2'    ))
  deo      = dataset.createVariable('deo'     , 'f8', ('deo_dim1'     , 'deo_dim2'    ))
  inlimits = dataset.createVariable('inlimits', 'f8', ('inlimits_dim1','inlimits_dim2'))
      
  weights [:] = Tweights 
  conec   [:] = Tconec   
  units   [:] = Tunits   
  inno    [:] = Tinno    
  outno   [:] = Toutno   
  eni     [:] = Teni     
  deo     [:] = Tdeo     
  inlimits[:] = Tinlimits

  dataset.close()
  
  return


def Print_Header ( LUT_file=None, coeff_type=None, ANN_basename=None, connectivity=None, test_perc=None, err_inc=None, 
                   nproc=None, shuffle_alw=None, netcdf=None, info=None ):
  print '################# Calculating ANN #################'
  print '                                                   '
  print 'LUT_file      ::  {}'.format(LUT_file)
  print 'coeff_type    ::  {}'.format(coeff_type)
  print 'ANN_basename  ::  {}'.format(ANN_basename)
  print 'connectivity  ::  {}'.format(connectivity)
  print 'test_perc     ::  {}'.format(test_perc)
  print 'err_inc       ::  {}'.format(err_inc)
  print 'nproc         ::  {}'.format(nproc)
  print 'shuffle_alw   ::  {}'.format(shuffle_alw)
  print 'netcdf        ::  {}'.format(netcdf)
  print '                                                   '
  if info!=None:
    print '---------------------------------------------------'
    print '                                                   '
    print '{}'.format(info)
    print '                                                   '
  print '###################################################'
  return


def Calc_ANN ( LUT_file, coeff_type, ANN_setup, ANN_basename='network', connectivity='standard', test_perc=0.3, err_inc=0.01,
               nproc=8, shuffle_alw=True, info=None ):
  from ffnet import savenet, loadnet

  if type(ANN_setup)==str:
    net = loadnet(ANN_setup)
    connectivity='UNKNOWN'; ANN_setup=ANN_setup+' is used ...'
  elif type(ANN_setup)==tuple or type(ANN_setup)==list:
    net = Init_ANN ( ANN_setup, connectivity=connectivity )
  else:
    raise ValueError('unknow ANN_setup input; use an existing ANN (saved with ffnet.savenet) or the structure for a new one')
  
  Print_Header ( LUT_file=LUT_file, coeff_type=coeff_type, ANN_basename=ANN_basename, connectivity=connectivity, 
                 test_perc=test_perc, err_inc=err_inc, nproc=nproc, shuffle_alw=shuffle_alw,
                 info=info )

  src, trg = Getting_Arrays ( LUT_file, coeff_type )
  indices, s_src = Shuffle_2D_X ( src )
  s_trg = trg[indices,:]; del(indices)
  
  net, err, rel_err = Train_ANN ( net, s_src, s_trg, test_perc=test_perc, nproc=nproc )
  print 'err:\terr_train:\terr_test:\trel_err:\trel_err_train:\trel_err_test:'
  print '{}\t{}\t{}\t{}\t{}\t{}'.format(err[0],err[1],err[2],rel_err[0],rel_err[1],rel_err[2])
  savenet(net, ANN_basename + '_0_.net')

  if type(LUT_file)==list:
    sl = LUT_file[0].rfind('/')
    if sl<0:
      path = ''
    else:
      path = LUT_file[0][:sl]
  else:
    sl = LUT_file.rfind('/')
    if sl<0:
      path = ''
    else:
      path = LUT_file[:sl]

  if coeff_type=='diff2diff' or coeff_type='diffuse':
    netcdf = path+'LUT_dstorder_8_10.diffuse.dx67.pspace.dz20.kabs20.ksca20.g3.delta_T_1.000_diff2diff.ANN.nc'
  elif coeff_type=='dir2diff':
    netcdf = path + 'LUT_dstorder_8_10.direct.dx67.pspace.dz20.kabs20.ksca20.g3.phi10.theta19.delta_T_1.000_dir2diff.ANN.nc'
  elif coeff_type=='dir2dir':
    netcdf = path + 'LUT_dstorder_8_10.direct.dx67.pspace.dz20.kabs20.ksca20.g3.phi10.theta19.delta_T_1.000_dir2dir.ANN.nc'
  else:
    raise ValueError ( 'bled' )

  ANN_to_NetCDF(net, netcdf, iprint=False)
  
  num, err_change = 1, -1.0
  while err_change<=err_inc:

    if shuffle_alw:
      indices, s_src = Shuffle_2D_X ( src )
      s_trg = trg[indices,:]
    
    err_old = rel_err[2]/100.0
    net, err, rel_err = Train_ANN ( net, s_src, s_trg, test_perc=test_perc, nproc=nproc )
    err_change = (rel_err[2]/100.0)-err_old
    print '{}\t{}\t{}\t{}\t{}\t{}'.format(err[0],err[1],err[2],rel_err[0],rel_err[1],rel_err[2])
    
    savenet(net, ANN_basename+'_{}_.net'.format(num)); num=num+1
    ANN_to_NetCDF(net, netcdf, iprint=False)

  print '\n{}-ANN was successfully calculated!!!'.format(coeff_type)
  print 'final err_change :: {}'                 .format(err_change)
  
  return
