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


def Getting_Arrays_Diffuse ( LUT_file ):
  
  import netCDF4 as nc
  import numpy as np

  try:
    ind = LUT_file.rfind('dx')+2; i=0 
    while LUT_file[ind+i]!='.': i=i+1
    dx = int(LUT_file[ind:ind+i])

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

  except:
    print 'Error in "Getting_Arrays_Diffuse"; probably error when reading netCDF data'
    print 'Following variables exist in the LUT-file:'
    print LUT_data.variables.keys()

  return index, src, S, S_tol


def Get_Dim ( LUT_file, key ):
  
  ind = LUT_file.rfind( key ) + len(key); num=0
  while LUT_file[ind+num]!='.': num=num+1
  dim = int(LUT_file[ind:ind+num])

  return dim


def Getting_Arrays_Direct ( LUT_file ):

  import netCDF4 as nc
  import numpy as np

# try:
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

# except:
#   print 'Error in "Getting_Arrays_Direct"; probably error when reading netCDF data'
#   print 'Following variables exist in the LUT-file:'
#   print LUT_data.variables.keys()
  
  return index, src, S, S_tol, T, T_tol


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

  try:
    net.train_tnc ( src_train, trg_train, nproc='ncpu' )
  except:
    net.train_tnc ( src_train, trg_train, nproc=nproc  )

  name = ['complete','train','test']; err, rel_err = [], []
  for i in range(len(name)):
    out, _ = net.test ( src_list[i], trg_list[i], iprint=0 )
    xerr, xrel_err = RMSE ( out, trg_list[i] )
    err.append([xerr, name[i]]); rel_err.append([xrel_err, name[i]])
  
  return net, err, rel_err


def ANN_to_NetCDF ( net_file, out_file='neural_network.nc' ):
    
  import netCDF4 as nc
  import numpy as np
  import ffnet as ff    
  
  print 'Exporting Network "%s" to .nc file :: %s' %(net_file, out_file)

  network = ff.loadnet ( net_file )
  
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


def Calc_ANN_Diffuse ( LUT_file, ANN_setup, connectivity='standard', test_perc=0.3, err_inc=1.0, nproc=8,
                       shuffle_alw=True, ANN_basename=None, netcdf_ANN=None ):
  
  from ffnet import savenet

  _, src, S, _ = Getting_Arrays_Diffuse ( LUT_file )

  indices, s_src = Shuffle_2D_X ( src )
  s_S = S[indices,:]
  
  print '##########  Calculating Diffuse ANN   ############'
  print '                                                  '
  print 'setup     ::    {}          '.format(ANN_setup)
  print 'ANN_type  ::    {}-connected'.format(connectivity)
  print 'test_perc ::    {}          '.format(test_perc)
  print 'err_inc   ::    {} (in %)   '.format(err_inc)
  print 'nproc     ::    {}/"ncpu"   '.format(nproc)
  print 'LUT       ::    {}          '.format(LUT_file)
  print '                                                  '
  print '##################################################'
  print '                                                  '


  net = Init_ANN ( ANN_setup, connectivity=connectivity )

  net, err, rel_err = Train_ANN ( net, s_src, s_S, test_perc=test_perc, nproc=nproc )

  print 'err:\terr_train:\terr_test:\trel_err:\trel_err_train:\trel_err_test:'
  print '{}\t{}\t{}\t{}\t{}\t{}'.format(err[0][0],err[1][0],err[2][0],rel_err[0][0],rel_err[1][0],rel_err[2][0])
  if ANN_basename!=None:
    savenet(net, ANN_basename + '_0_.net')
  
  num=1
  while err_change<=err_inc:

    if shuffle_alw:
      indices, s_src = Shuffle_2D_X ( src )
      s_S[:,:] = S[indices,:]

    err_old = rel_err[2][0]*1.0
    net, err, rel_err = Train_ANN ( net, s_src, s_S, test_perc=test_perc, nproc=nproc )

    err_change = rel_err[2][0]-err_old

    print '{}\t{}\t{}\t{}\t{}\t{}'.format(err[0][0],err[1][0],err[2][0],rel_err[0][0],rel_err[1][0],rel_err[2][0])
    if ANN_basename!=None: ff.savenet(net, ANN_basename + '_{}_.net'.format(num))
    num=num+1

  if netcdf_ANN!=None:
    ANN_to_NetCDF ( net, netcdf_ANN )
    print 'ANN converted to netCDF format :: {}'.format(netcdf_ANN)
  
  print 'Diffuse ANN successfully calculated!!!'
  print 'Final err_change :: {}'.format(err_change)
  
  return net


def Calc_ANN_Direct ( LUT_file, ANN_setup, LUT_type, connectivity='standard', test_perc=0.3, err_inc=1.0, 
                      nproc=8, shuffle_alw=True, ANN_basename=None, netcdf_ANN=None ):

  from ffnet import savenet
  from numpy import append

  _, src, S, _, T, _ = Getting_Arrays_Direct ( LUT_file[0] )
  if len(LUT_file)>1:
    for LUT in range(len(LUT_file)-1):
      _, xsrc, xS, _, xT, _ = Getting_Arrays_Direct ( LUT_file[LUT+1] )
      src = append(src, xsrc, axis=0)
      S   = append(S  , xS  , axis=0)
      T   = append(T  , xT  , axis=0)

  indices, s_src = Shuffle_2D_X ( src )

  if LUT_type=='dir2dir' or LUT_type=='T':
    s_trg = T[indices,:]; trg=T*1.0; del(S); del(T)
  elif LUT_type=='dir2diff' or LUT_type=='S':
    s_trg = S[indices,:]; trg=S*1.0; del(S); del(T)
  else:
    raise ValueError ('unknown "LUT_type" :: {}\t only "dir2dir"/"T" or "dir2diff"/"S" are possible'.format(LUT_type))

  print '###########  Calculating Direct ANN   ############'
  print '                {}-network  '.format(LUT_type)
  print '                            '
  print 'setup     ::    {}          '.format(ANN_setup)
  print 'ANN_type  ::    {}-connected'.format(connectivity)
  print 'test_perc ::    {}          '.format(test_perc)
  print 'err_inc   ::    {} (in %)   '.format(err_inc)
  print 'nproc     ::    {}/"ncpu"   '.format(nproc)
  print 'LUT       ::    {}          '.format(LUT_file)
  print '                                                  '
  print '##################################################'
  print '                                                  '

  net = Init_ANN ( ANN_setup, connectivity=connectivity )

  net, err, rel_err = Train_ANN ( net, s_src, s_trg, test_perc=test_perc, nproc=nproc )

  print 'err:\terr_train:\terr_test:\trel_err:\trel_err_train:\trel_err_test:'
  print '{}\t{}\t{}\t{}\t{}\t{}'.format(err[0][0],err[1][0],err[2][0],rel_err[0][0],rel_err[1][0],rel_err[2][0])
  if ANN_basename!=None:
    savenet(net, ANN_basename + '_0_.net')
  
  num=1; err_change=-1.0
  while err_change<=err_inc:

    if shuffle_alw:
      indices, s_src = Shuffle_2D_X ( src )
      s_trg[:,:] = trg[indices,:]

    err_old = rel_err[2][0]*1.0
    net, err, rel_err = Train_ANN ( net, s_src, s_trg, test_perc=test_perc, nproc=nproc )

    err_change = rel_err[2][0]-err_old

    print '{}\t{}\t{}\t{}\t{}\t{}'.format(err[0][0],err[1][0],err[2][0],rel_err[0][0],rel_err[1][0],rel_err[2][0])
    if ANN_basename!=None: savenet(net, ANN_basename + '_{}_.net'.format(num))
    num=num+1

  if netcdf_ANN!=None:
    ANN_to_NetCDF ( net, netcdf_ANN )
    print 'ANN converted to netCDF format :: {}'.format(netcdf_ANN)
  
  print 'Direct ({}) ANN successfully calculated!!!'.format(LUT_type)
  print 'Final err_change :: {}'.format(err_change)

  return net
