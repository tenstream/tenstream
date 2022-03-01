#!/usr/bin/env python
# -*- coding: utf-8 -*-

import netCDF4 as nc
import numpy as np
import ffnet as ff
import argparse


parser = argparse.ArgumentParser(prog='Program which calculates and trains Artificial Neural Networks (ANN) to fit tenstream-Look-Up-Tables (LUT)')
parser.add_argument('-l', '--LUT', nargs='+', type=str, required=True, help='LUT/LUTs which are used to train')
parser.add_argument('-s', '--ANN_setup', nargs='+', type=str, required=True, 
                    help='number of hidden neurons and hidden layers (e.g. use 20 20 for a (input_nodes,20,20,output_nodes)-network) OR existing ANN which has been saved with "ffnet.savenet"')
parser.add_argument('-c', '--coeff_type', choices=['diffuse','diff2diff','dir2diff','dir2dir'], type=str, required=True, 
                    help='use diff2diff, dir2diff, or dir2dir; if you want to calculate a diff2diff/dir2* network you have to use a diffuse/direct LUT for training')
parser.add_argument('-t', '--test_perc', type=float, required=True, help='percentage of training data which is used to test')
parser.add_argument('-i', '--err_inc', type=float, required=True, help='maximal error increase of test data; a larger error increase leads to a stop of the training')
parser.add_argument('-b', '--basename', type=str, required=True, help='ffnet is saved to basename_<training step>_.net after every training step')
parser.add_argument('--full', action='store_true', help='if set a fully connected ANN will be initialized and trained')
parser.add_argument('--nproc', default=8, help='number of processes; use ncpu for maximal possible number')
parser.add_argument('--index', action='store_true', help='calc an "index-ANN"')
args = parser.parse_args()


# process input:
#----------------------------------------------------
if args.coeff_type=='dir2dir':
    LUT_var = 'T' 
else:
    LUT_var = 'S'


try:
    nproc = int(args.nproc)
except ValueError:
    nproc = args.nproc

try:
    args_ANN_setup = map(int, args.ANN_setup)
except ValueError:
    args_ANN_setup = args.ANN_setup[0]
#----------------------------------------------------


#----------------------------------------------------
################## functions ########################
#----------------------------------------------------

def RMSE( a, b, weights=None ):
    
    if np.shape(a)!=np.shape(b): raise ValueError('Arrays do not have same shapes :: {} .ne. {}'.format(np.shape(a),np.shape(b)) )
    s= np.sqrt ( np.average( (a - b)**2, weights=weights ) )
      
    return np.array([s,s/np.maximum(1e-8,abs(np.average(b,weights=weights)))*100])


def Shuffle_2D_X ( *arr ):
    
    indices = np.arange(0, len(arr[0]), dtype=int)
    np.random.shuffle(indices)
    s_arr = [xarr[indices,:] for xarr in arr]
    if len(s_arr)==1: s_arr=s_arr[0]
    
    return s_arr


def Get_Dim ( LUT_file, *key ):
   
    dim = {}
    for xkey in key:
        ind = LUT_file.rfind( xkey ) + len(xkey); num=0
        while LUT_file[ind+num]!='.': num=num+1
        dim[xkey] = int(LUT_file[ind:ind+num])

    return dim


def Getting_Arrays ( LUT_file ):
    
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
        var.extend( ['phi','theta'] )

        phi   = LUT_data.variables['direct.dx{0:d}.dy{1:d}.pspace.phi'  .format(res['dx'],res['dy'])][:].tolist()
        theta = LUT_data.variables['direct.dx{0:d}.dy{1:d}.pspace.theta'.format(res['dx'],res['dy'])][:].tolist()
    

    data = [LUT_data.variables['{0:s}.dx{1:d}.dy{2:d}.pspace.{3:s}'.format(coeff_type,res['dx'],res['dy'],xvar)][:] for xvar in var]
    LUT  = [LUT_data.variables['{0:s}.{1:s}.{2:s}'.format(coeff_type,LUT_prefix,LUT_name)][:] for LUT_name in LUT_names]
    
    for  dz  in range( res[ 'dz' ] ):
        for kabs in range( res['kabs'] ):
            for ksca in range( res['ksca'] ):
                for  g   in range( res[ 'g'  ] ):

                    src  .append ( [data[0][dz], data[1][kabs], data[2][ksca], data[3][g]] )
                    S    .append ( LUT[0][g,ksca,kabs,dz] )
                    S_tol.append ( LUT[1][g,ksca,kabs,dz] )
                    if coeff_type=='direct':
                        T    .append ( LUT[2][g,ksca,kabs,dz] )
                        T_tol.append ( LUT[3][g,ksca,kabs,dz] )
                        index.append ( [dz,kabs,ksca,g,phi.index(float(res['phi'])),theta.index(float(res['theta']))] )
                    else:
                        index.append ( [dz,kabs,ksca,g] )

    dic = {'index':np.array(index), 'S':np.array(S), 'S_tol':np.array(S_tol)}; src=np.array(src)
    for ind, xvar in enumerate(var): dic[xvar]=np.array(data[ind])
    if coeff_type=='direct': 
        dic.update( {'T':np.array(T), 'T_tol':np.array(T_tol)} )
        app = np.append ( np.ones((len(src),1),dtype=float)*res['phi'], np.ones((len(src),1),dtype=float)*res['theta'], axis=1 )
        src = np.append ( src, app, axis=1 )

    dic.update( {'src':src} )


    return dic


def Get_Output_Name ( LUT_fname, coeff_type ):

    path     = LUT_fname[:LUT_fname.rfind('LUT')]
    fname    = LUT_fname[LUT_fname.rfind('LUT'):]
    dim      = Get_Dim ( LUT_fname, 'dx', 'dz', 'kabs', 'ksca', 'g' )
    out_name = path+'LUT_dstorder_8_10.diffuse.dx{0:d}.pspace.dz{1:d}.kabs{2:d}.ksca{3:d}.g{4:d}.' \
                    .format(dim['dx'],dim['dz'],dim['kabs'],dim['ksca'],dim['g']) 
    if fname.rfind('direct')>=0: 
        dim.update( {'phi':10, 'theta':19} )                   # hardcoded
        out_name = out_name.replace('diffuse','direct')
        out_name = out_name+'phi{0:d}.theta{1:d}.'.format(dim['phi'],dim['theta'])
    
    return out_name+'delta_T_1.000_{0:s}.ANN.nc'.format(coeff_type)


def ANN_to_NetCDF ( net, out_file, iprint=True, **data ):
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
    for key, val in data.iteritems(): data[key] = val.T

    dataset = nc.Dataset(out_file, 'w', format='NETCDF4')
    
    dataset.createDimension ( 'weights_dim1'  , np.shape(Tweights )[0] )
    dataset.createDimension ( 'conec_dim1'    , np.shape(Tconec   )[0] ); dataset.createDimension ( 'conec_dim2'   , np.shape(Tconec   )[1] )
    dataset.createDimension ( 'units_dim1'    , np.shape(Tunits   )[0] )
    dataset.createDimension ( 'inno_dim1'     , np.shape(Tinno    )[0] )
    dataset.createDimension ( 'outno_dim1'    , np.shape(Toutno   )[0] )
    dataset.createDimension ( 'eni_dim1'      , np.shape(Teni     )[0] ); dataset.createDimension ( 'eni_dim2'     , np.shape(Teni     )[1] )
    dataset.createDimension ( 'deo_dim1'      , np.shape(Tdeo     )[0] ); dataset.createDimension ( 'deo_dim2'     , np.shape(Tdeo     )[1] )
    dataset.createDimension ( 'inlimits_dim1' , np.shape(Tinlimits)[0] ); dataset.createDimension ( 'inlimits_dim2', np.shape(Tinlimits)[1] )
    for key, val in data.iteritems(): 
        dataset.createDimension ( 'pspace.{}_dim1'.format(key), np.shape(data[key])[0])
  
    weights  = dataset.createVariable('weights' , 'f8',  'weights_dim1'                  )
    conec    = dataset.createVariable('conec'   , 'i' , ('conec_dim1'   , 'conec_dim2'  ))
    units    = dataset.createVariable('units'   , 'f8',  'units_dim1'                    )
    inno     = dataset.createVariable('inno'    , 'i' ,  'inno_dim1'                     )
    outno    = dataset.createVariable('outno'   , 'i' ,  'outno_dim1'                    )
    eni      = dataset.createVariable('eni'     , 'f8', ('eni_dim1'     , 'eni_dim2'    ))
    deo      = dataset.createVariable('deo'     , 'f8', ('deo_dim1'     , 'deo_dim2'    ))
    inlimits = dataset.createVariable('inlimits', 'f8', ('inlimits_dim1','inlimits_dim2'))
    dataset_list = {}
    for key, val in data.iteritems(): 
        dataset_list [key] = dataset.createVariable('pspace.{}'.format(key), 'f8', 'pspace.{}_dim1'.format(key))

    weights [:] = Tweights 
    conec   [:] = Tconec   
    units   [:] = Tunits   
    inno    [:] = Tinno    
    outno   [:] = Toutno   
    eni     [:] = Teni     
    deo     [:] = Tdeo     
    inlimits[:] = Tinlimits
    for key, var in dataset_list.iteritems(): var [:] = data[key]

    dataset.close()
    
    return


def Test_ANN ( net, *arr ):
    
    err, rel_err = [], []
    for xarr in arr:
        out, _ = net.test( xarr[0], xarr[1], iprint=0 )
        xerr, xrel_err = RMSE( out, xarr[1] )
        err.append(xerr); rel_err.append(xrel_err)
    err.extend(rel_err)

    return err


def Print_Header ( **info ):
    print '##################### ANN #####################'
    print '                                               '
    for key, value in info.iteritems():
        print '{}\t::\t{}'.format(key,value)
    print '                                               '
    print '###############################################'
    print '                                               '
    print '-----------------------------------------------'
    print '                                               '
    return
#----------------------------------------------------


#----------------------------------------------------
#################### main ###########################
#----------------------------------------------------

# load source and target array    
LUT = [Getting_Arrays(LUT_file) for LUT_file in args.LUT]
if args.index:
    src = np.concatenate([xLUT[ 'index' ] for xLUT in LUT], axis=0)
    src = src.astype( float )
else:
    src = np.concatenate([xLUT[ 'src'   ] for xLUT in LUT], axis=0)
trg = np.concatenate([xLUT[LUT_var] for xLUT in LUT], axis=0)


# initialize ANN
if type(args_ANN_setup)==str:
    net = ff.loadnet( args_ANN_setup )
    setup = "UNKNOWN"
else:
    num_inp_nodes, num_out_nodes = [src.shape[1]], [trg.shape[1]]
    setup = tuple(num_inp_nodes+args_ANN_setup+num_out_nodes)
    if args.full:
        conec = ff.tmlgraph( setup )
    else:
        conec = ff.mlgraph ( setup )
    net = ff.ffnet( conec )

# build up a shuffled test and train array
ind = int( src.shape[0]*(1.0-args.test_perc) )
s_src, s_trg = Shuffle_2D_X ( src, trg )
src_train, trg_train = s_src[:ind,:], s_trg[:ind,:]
src_test , trg_test  = s_src[ind:,:], s_trg[ind:,:]

# output
netcdf_output = Get_Output_Name( args.LUT[0], args.coeff_type )
Print_Header ( LUT_name=args.LUT, coeff_type=args.coeff_type, ANN_setup=setup, test_percentage=args.test_perc, num_of_proc=args.nproc,
               basename=args.basename, netcdf_output=netcdf_output, err_increase=args.err_inc )

# train the network for the first time
net.train_tnc ( src_train, trg_train, nproc=nproc )
err = Test_ANN ( net, [src,trg], [src_train,trg_train], [src_test,trg_test] )
print 'err:\terr_train:\terr_test:\trel_err [%]:\trel_err_train [%]:\trel_err_test [%]:'
print '{}\t{}\t{}\t{}\t{}\t{}'.format(*err)

# save network as ffnet and in netCDF4-format:
ff.savenet( net, args.basename + '_1_.net' ); train_num = 2
if args.coeff_type=='diffuse' or args.coeff_type=='diff2diff':
    ANN_to_NetCDF ( net, netcdf_output, iprint=False, dz=LUT[0]['dz'], kabs=LUT[0]['kabs'], ksca=LUT[0]['ksca'], g=LUT[0]['g'] )
else:
    ANN_to_NetCDF ( net, netcdf_output, iprint=False, dz=LUT[0]['dz'], kabs=LUT[0]['kabs'], ksca=LUT[0]['ksca'], g=LUT[0]['g'],
                    phi=LUT[0]['phi'], theta=LUT[0]['theta'] )

# train the network until the error of the test data increases more than "err_inc"
err_test_ch = -999.9
while err_test_ch<args.err_inc:
    err_old_test = err[-1]/100.0

    s_src, s_trg = Shuffle_2D_X ( src, trg )
    src_train, trg_train = s_src[:ind,:], s_trg[:ind,:]
    src_test , trg_test  = s_src[ind:,:], s_trg[ind:,:]

    net.train_tnc ( src_train, trg_train, nproc=nproc )
    
    err = Test_ANN ( net, [src, trg], [src_train, trg_train], [src_test, trg_test] )
    err_test_ch  = err[-1]/100.0 - err_old_test
    print '{}\t{}\t{}\t{}\t{}\t{}'.format(*err)

    ff.savenet( net, args.basename + '_{}_.net'.format(train_num) ); train_num=train_num+1
    if args.coeff_type=='diffuse' or args.coeff_type=='diff2diff':
        ANN_to_NetCDF ( net, netcdf_output, iprint=False, dz=LUT[0]['dz'], kabs=LUT[0]['kabs'], ksca=LUT[0]['ksca'], g=LUT[0]['g'] )
    else:
        ANN_to_NetCDF ( net, netcdf_output, iprint=False, dz=LUT[0]['dz'], kabs=LUT[0]['kabs'], ksca=LUT[0]['ksca'], g=LUT[0]['g'],
                        phi=LUT[0]['phi'], theta=LUT[0]['theta'] )
