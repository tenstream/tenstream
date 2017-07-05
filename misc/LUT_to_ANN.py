#!/usr/bin/env python

import numpy as np
import netCDF4 as NC
import sys
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm

def rmse(a,b,weights=None):
        import numpy as np
        if np.shape(a)!=np.shape(b): raise ValueError('Arrays do not have same shapes :: {} .ne. {}'.format(np.shape(a),np.shape(b)) )

        s= np.sqrt ( np.average( (a - b)**2, weights=weights ) )
        return np.array([s,s/np.maximum(1e-8,abs(np.average(b,weights=weights)))*100])

def export_ANN_to_nc ( network, filename, Nlayers, Nneurons, pspace):

    import netCDF4 as nc
    import numpy as np
    import os as os

    print('Exporting Network to .nc file :: ',filename)

    if os.path.isfile(filename):
        os.remove(filename)

    # transpose arrays to fortran order
    Tweights  = network.weights .T
    Tconec    = network.conec   .T
    Tunits    = network.units   .T
    Tinno     = network.inno    .T
    Toutno    = network.outno   .T
    Teni      = network.eni     .T
    Tdeo      = network.deo     .T
    Tinlimits = network.inlimits.T
    for key, val in pspace.items():
      pspace[key] = val.T

    try:
        dataset = nc.Dataset(filename, 'w')

        dataset.createDimension ( 'weights_dim1' , np.shape(Tweights )[0] )
        dataset.createDimension ( 'conec_dim1'   , np.shape(Tconec   )[0] ); dataset.createDimension ( 'conec_dim2'   , np.shape(Tconec   )[1] )
        dataset.createDimension ( 'units_dim1'   , np.shape(Tunits   )[0] )
        dataset.createDimension ( 'inno_dim1'    , np.shape(Tinno    )[0] )
        dataset.createDimension ( 'outno_dim1'   , np.shape(Toutno   )[0] )
        dataset.createDimension ( 'eni_dim1'     , np.shape(Teni     )[0] ); dataset.createDimension ( 'eni_dim2'     , np.shape(Teni     )[1] )
        dataset.createDimension ( 'deo_dim1'     , np.shape(Tdeo     )[0] ); dataset.createDimension ( 'deo_dim2'     , np.shape(Tdeo     )[1] )
        dataset.createDimension ( 'inlimits_dim1', np.shape(Tinlimits)[0] ); dataset.createDimension ( 'inlimits_dim2', np.shape(Tinlimits)[1] )

        for key, val in pspace.items():
            dataset.createDimension ( 'pspace.{}_dim1'.format(key), np.shape(pspace[key])[0])


        weights  = dataset.createVariable('weights' , 'f8',  'weights_dim1'                  )
        conec    = dataset.createVariable('conec'   , 'i' , ('conec_dim1'   , 'conec_dim2'  ))
        units    = dataset.createVariable('units'   , 'f8',  'units_dim1'                    )
        inno     = dataset.createVariable('inno'    , 'i' ,  'inno_dim1'                     )
        outno    = dataset.createVariable('outno'   , 'i' ,  'outno_dim1'                    )
        eni      = dataset.createVariable('eni'     , 'f8', ('eni_dim1'     , 'eni_dim2'    ))
        deo      = dataset.createVariable('deo'     , 'f8', ('deo_dim1'     , 'deo_dim2'    ))
        inlimits = dataset.createVariable('inlimits', 'f8', ('inlimits_dim1','inlimits_dim2'))

        dataset_list = {}
        for key, val in pspace.items():
            dataset_list [key] = dataset.createVariable('pspace.{}'.format(key), 'f8', 'pspace.{}_dim1'.format(key))


        weights [:] = Tweights
        conec   [:] = Tconec
        units   [:] = Tunits
        inno    [:] = Tinno
        outno   [:] = Toutno
        eni     [:] = Teni
        deo     [:] = Tdeo
        inlimits[:] = Tinlimits

        for key, var in dataset_list.items():
            var [:] = pspace[key]

        dataset.nrhiddennodes = len(network.hidno)
        dataset.Nlayers  = Nlayers
        dataset.Nneurons = Nneurons

        dataset.close()
    except Exception as e:
        print('Error occured when we tried to export Network to file :: ',e)
        raise(e)

def scatter_to_hist(x,y,nbins=1000):
  import numpy as np
  hist,xedges,yedges = np.histogram2d(x.flatten(),y.flatten(),bins=nbins)
  extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
  return hist.T,extent

def plot_coeffs(ANNname, net,input,target,output):
    from pylab import figure,subplot,plot,title,tight_layout,savefig,matplotlib,imshow,colorbar
    import os
    print('plotting coefficients')
    fig = figure(figsize=(35,25) )
    Nstreams = np.minimum(149,np.shape(target)[1])
    nrows = int(np.sqrt(Nstreams)) ; ncols = int(np.ceil(1.*Nstreams/nrows))

    output = np.maximum(0., np.minimum(1., output ) )

    for i in np.arange(Nstreams):
            subplot(nrows,ncols,i+1)
            print('.')

            h,extent = scatter_to_hist(target[:,i],output[:,i],nbins=50)
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list('white to red over green', [(0,'white'),  (.01,'yellow'),  (.1,'green'), (.3,'orange'),  (.5,'red'), (1.,'purple')])
            imshow(h,origin='lower',extent=extent,cmap=cmap,norm=LogNorm())
            colorbar()

            v=np.linspace(np.min([target[:,i],output[:,i]]) , np.max( [target[:,i],output[:,i]] ) )
            plot(v,v,color='black',linestyle='-',linewidth=1, rasterized=True)
            ubound = np.minimum( v+5e-3, v*(1+.3) )
            lbound = np.maximum( v-5e-3, v*(1-.3) )
            plot(v,ubound,color='orange',linestyle='--',linewidth=1)
            plot(v,lbound,color='orange',linestyle='--',linewidth=1)
            error = np.around( rmse(target[:,i],output[:,i])[1] ,decimals=1)
            title(str(i)+'::'+str(error ))
            tight_layout()
    savefig('/tmp/network_plot_{}.pdf'.format(os.path.basename(ANNname)))
    return


def create_net(inno,layer_definitions,outno,network_function):
    from ffnet import ffnet
    pts = [inno,] + layer_definitions + [outno,]
    conec = network_function(pts, biases=True) ; net = ffnet(conec)
    print('created new network with ',pts,'layers and', len(net.conec),'connections')
    return net

def import_network( ANNname, Nlayers, Nneurons ):
    from ffnet import loadnet
    import os
    try:
        f = '{0:}_{1:}_{2:}'.format( ANNname,Nlayers,Nneurons)
        net = loadnet(f)
        print('Loaded network from file: ',f)
        print('Network limits:',net.inlimits)
        return net

    except Exception as e:
        print('Error when loading ffnet network file',e)
        raise(e)

def export_network( net, ANNname, Nlayers, Nneurons, pspace):
    from ffnet import savenet
    import os
    try:
        f = '{0:}_{1:}_{2:}'.format( ANNname, Nlayers, Nneurons)
        savenet(net, f)
        export_ANN_to_nc( net, ANNname+'.nc', Nlayers, Nneurons, pspace)

    except Exception as e:
        print('Error when writing ffnet network file',e)
        raise(e)


def compare_to_old_net(ANNname, Nlayers, Nneurons, new_net, test_inp, test_target, pspace):
  try:
    old_net = import_network(ANNname, Nlayers, Nneurons)
  except IOError as e:
    rmse_old_network = np.Inf
  else:
    output_old_network,_ = old_net.test(test_inp,test_target,iprint=0)
    rmse_old_network = rmse(test_target,output_old_network)[1]

  test_output,regression = new_net.test(test_inp,test_target,iprint=0)
  rmse_new_network = rmse(test_target,test_output)[1]

  print('Old network error: {0:} :: new network: {1:}'.format(rmse_old_network,rmse_new_network))

  if rmse_new_network < rmse_old_network:
      try:
          print('Saving new network to :::  ',ANNname)
          export_network( new_net, ANNname, Nlayers, Nneurons, pspace)

      except Exception as e:
          print('Error occured when we tried to save the ANN :: ',e)
          raise(e)
  else:
      print('Already existing, old network has already lower error.. old: {} new: {}'.format(rmse_old_network,rmse_new_network))

def train_net(ANNname, Nlayers, Nneurons, net, train_inp, train_target, test_inp, test_target, pspace, ncpu='ncpu'):
  print('training net: nr. of hidno: {} nr. of weights: {} using {} input entries'.format( len(net.hidno), np.shape(net.weights),np.shape(train_inp) ))

  last_weights    = net.weights
  last_train_err  = np.Inf
  last_test_err   = np.Inf
  iter            = 0
  index_successfull_update = 0

#TNC RC Strings  INFEASIBLE   = -1 # Infeasible (low > up)
#TNC RC Strings  LOCALMINIMUM =  0 # Local minima reach (|pg| ~= 0)
#TNC RC Strings  FCONVERGED   =  1 # Converged (|f_n-f_(n-1)| ~= 0)
#TNC RC Strings  XCONVERGED   =  2 # Converged (|x_n-x_(n-1)| ~= 0)
#TNC RC Strings  MAXFUN       =  3 # Max. number of function evaluations reach
#TNC RC Strings  LSFAIL       =  4 # Linear search failed
#TNC RC Strings  CONSTANT     =  5 # All lower bounds are equal to the upper bounds
#TNC RC Strings  NOPROGRESS   =  6 # Unable to progress
#TNC RC Strings  USERABORT    =  7 # User requested end of minimization

  while True:
    print('\n\n')
    maxiter = int( np.minimum( 1e3, np.maximum(1e2, len(net.weights)*2) ) )
#    import ipdb
#    ipdb.set_trace()

    try:
        res = net.train_tnc(train_inp,train_target,nproc=ncpu,maxfun=maxiter,disp=5)
#        res = net.train_basinhopping(train_inp,train_target,niter=maxiter,disp=5)
#        res = net.train_rprop(train_inp,train_target,maxiter=maxiter,disp=99)
#        res = net.train_bfgs(train_inp,train_target,bounds=( (None,None), )*len(net.weights),disp=0,maxfun=maxiter)
#        res = net.train_momentum(train_inp,train_target,disp=1,maxiter=100*len(net.weights))
#        res = net.train_genetic(train_inp,train_target,verbosity=0)
    except Exception as e:
        print('Error occured in ffnet training:',e)
        import traceback
        print(traceback.format_exc())

        print('\n\n ', type(train_inp), type(train_target))
        raise(e)

    train_output,regression = net.test(train_inp,train_target,iprint=0)

    test_output,regression = net.test(test_inp,test_target,iprint=0)

    iter = iter+1
    print(len(net.hidno),'::',iter,'RMSE train',rmse(train_target,train_output),' test ',rmse(test_target,test_output),'last:',[last_train_err,last_test_err])

    lbreak = False
    if rmse(test_target,test_output)[1] > 1.001*last_test_err:
      lbreak = True
      print('error for testing set, grew by',rmse(test_target,test_output)[1]/last_test_err * 100,'%')

    if iter - index_successfull_update > 20:
      lbreak=True
      print('tried to enhance result too often, its probably no use....')

    if lbreak:
      print('::','training ended!!!! -- stopping training!')
      compare_to_old_net(ANNname,Nlayers, Nneurons, net, test_inp, test_target, pspace)
      return

    if rmse(test_target,test_output)[1] < last_test_err:
      print('::','Updating last err. Err Increment: {}'.format(last_test_err-rmse(test_target,test_output)[1]))
      last_test_err = rmse(test_target,test_output)[1]
      last_weights = net.weights
      compare_to_old_net(ANNname,Nlayers, Nneurons, net, test_inp, test_target, pspace)
      index_successfull_update=iter

    last_train_err = rmse(train_target,train_output)[1]

def create_training_dataset(coeff_mode, LUTname, training_fraction=.8):
    # Try to find out which LUT entries we have here
    try:
        D = NC.Dataset(LUTname,'r')
        diffuse_tables = {k for k, v in D.variables.items() if k.endswith('.S') }
        direct_tables  = {k for k, v in D.variables.items() if k.endswith('.T') }

        for table in diffuse_tables: print(' Found diffuse Table:', table)
        for table in direct_tables:  print(' Found direct  Table:', table)

        for table in diffuse_tables:
            # try to find sza and azi angle from table name
            dx = dy = phi = theta = None
            descr = table.replace('/','.').split('.')
            for k in descr:
                if(k.find('dx')!=-1):    dx    = k[2:]
                if(k.find('dy')!=-1):    dy    = k[2:]
                if(k.find('phi')!=-1):   phi   = k[3:]
                if(k.find('theta')!=-1): theta = k[5:]

    except Exception as e:
        print('Error occured reading LUTfile: ',LUTname,' because',e)
        raise(e)
    finally:
        D.close()


    if coeff_mode=='diff2diff':
        if len(diffuse_tables)!=1:
            raise(Exception('I expect exactly one diffuse LUT in file -- what is this? I found {0:} diffuse tables for coefficients'.format(len(diffuse_tables))) )

    if coeff_mode=='dir2dir':
        if len(direct_tables)<1:
            raise(Exception('I expect at least one direct LUT in file -- what is this? I found {0:} direct tables for coefficients'.format(len(direct_tables))) )

    if coeff_mode=='dir2diff':
        if len(direct_tables)<1:
            raise(Exception('Although I dont use direct tables for dir2diff, I expect at least one direct LUT in the file -- something out of convention? I found {0:} direct tables for coefficients'.format(len(direct_tables))) )


    print('\n Creating training and testing datasets:')

    inp = []
    out = []
    pspace = {}

    if coeff_mode=='diff2diff':
        try:
            D = NC.Dataset(LUTname,'r')
            table = [k for k, v in D.variables.items() if k.endswith('.S') ][0]
            # try to find sza and azi angle from table name
            dx = dy = phi = theta = None
            descr = table.replace('/','.').replace('.S','')
            for k in descr.split('.'):
                if(k.find('dx')!=-1):    dx    = k[2:]
                if(k.find('dy')!=-1):    dy    = k[2:]
                if(k.find('phi')!=-1):   phi   = k[3:]
                if(k.find('theta')!=-1): theta = k[5:]

            S     = D.variables[table][:].T
            g     = D.variables[descr+'.pspace.g'][:]
            aspect= D.variables[descr+'.pspace.aspect'][:]
            tau   = D.variables[descr+'.pspace.tau'][:]
            w0    = D.variables[descr+'.pspace.w0'][:]
            pspace['g'     ] = g
            pspace['aspect'] = aspect
            pspace['tau'   ] = tau
            pspace['w0'    ] = w0

        except Exception as e:
            print('Error occured reading LUTfile: ',LUTname,' because',e)
            raise(e)
        finally:
            D.close()

        for iaspect in range(np.size(aspect)):
            for itau in range(np.size(tau)):
                for iw0 in range(np.size(w0)):
                    for ig in range(np.size(g)):
                        # inp.append( [ dz[idz],kabs[ikabs],ksca[iksca],g[ig] ] )
                        inp.append( [ iaspect, itau, iw0, ig ] )
                        out.append( S[:, iaspect, itau,iw0, ig] )

    if coeff_mode=='dir2dir' or coeff_mode=='dir2diff':
        try:
            D = NC.Dataset(LUTname,'r')

            # loop over all LUT entries in file:
            if coeff_mode=='dir2diff':
                tables = {k for k, v in D.variables.items() if k.endswith('.S') }

            if coeff_mode=='dir2dir':
                tables  = {k for k, v in D.variables.items() if k.endswith('.T') }

            for table in tables:
                # try to find sza and azi angle from table name
                dx = dy = phi = theta = None
                descr = table.replace('/','.').replace('.S','').replace('.T','')
                for k in descr.split('.'):
                    if(k.find('dx')!=-1):    dx    = k[2:]
                    if(k.find('dy')!=-1):    dy    = k[2:]
                    if(k.find('phi')!=-1):   phi   = k[3:]
                    if(k.find('theta')!=-1): theta = k[5:]

                print('Reading dataset from: {0:} with descr {1:}'.format(table,descr))

                basedescr = '.'.join(descr.split('.')[:1])

                SorT  = D.variables[table][:].T
                aspect= D.variables[basedescr+'.pspace.aspect'][:]
                tau   = D.variables[basedescr+'.pspace.tau'][:]
                w0    = D.variables[basedescr+'.pspace.w0'][:]
                g     = D.variables[basedescr+'.pspace.g'][:]

                pspace['g'     ] = aspect
                pspace['aspect'] = tau
                pspace['tau'   ] = w0
                pspace['w0'    ] = g
                pspace['phi'  ] = D.variables[basedescr+'.pspace.phi'][:]
                pspace['theta'] = D.variables[basedescr+'.pspace.theta'][:]

                for iaspect in range(np.size(aspect)):
                    for itau in range(np.size(tau)):
                        for iw0 in range(np.size(w0)):
                            for ig in range(np.size(g)):
                                # inp.append( [ dz[idz],kabs[ikabs],ksca[iksca],g[ig],phi,theta ] )
                                inp.append( [ iaspect, itau, iw0, ig, phi, theta ] )
                                out.append( SorT[:,iaspect,itau,iw0,ig] )

        except Exception as e:
            print('Error occured reading LUTfile: ',LUTname,' because',e)
            raise(e)
        finally:
            D.close()

    inp = np.array(inp); out = np.array(out)
    inp = inp.astype(np.float)
    out = out.astype(np.float)

    # Shuffle entries
    N_entries = np.shape(inp)[0]
    permu = np.random.permutation(N_entries)
    inp = inp[permu] ; out = out[permu]

    # up until entries to 'divider' are used for training, the rest for testing
    divider = int(np.round(N_entries*training_fraction))

    print(' using {} points for training and {} for testing'.format(divider, N_entries-divider))
    train_inp = inp[:divider,:] ; train_target = out[:divider,:]
    test_inp  = inp[divider:,:] ; test_target  = out[divider:,:]

    return np.array([train_inp, train_target, test_inp, test_target]), pspace

def train_ANN(ANNname, Nlayers, Nneurons, train_dataset, pspace, ncpu='ncpu'):
    from ffnet import mlgraph
    train_inp, train_target, test_inp, test_target = train_dataset
    try:
        net = import_network( ANNname, Nlayers, Nneurons)
    except Exception as e:
        print('Could not load existing neural network: {0:} -- creating a new one'.format(e))
        hiddenlayers = [ Nneurons ] * Nlayers
        net = create_net(np.shape(train_inp)[1], hiddenlayers, np.shape(train_target)[1], mlgraph)

    export_network( net, ANNname, Nlayers, Nneurons, pspace)

    compare_to_old_net(ANNname, Nlayers, Nneurons, net,test_inp,test_target, pspace)

    train_net(ANNname, Nlayers, Nneurons, net, train_inp, train_target, test_inp, test_target, pspace, ncpu=ncpu)

    output,_ = net.test(test_inp, test_target, iprint=0)
    plot_coeffs(ANNname, net, test_inp, test_target, output)



if __name__ == "__main__":
  import argparse
  import glob

  parser = argparse.ArgumentParser(description='Input parameters for Neural Network Generator')
  parser.add_argument('LUTfiles', type=str, nargs='+', help='LUT basename which should be converted')
  parser.add_argument('coeffmode', choices=['diff2diff', 'dir2dir', 'dir2diff'])
  parser.add_argument('ANNname', type=str, help='ANN file name on which the net will be saved')
  parser.add_argument('-N', '--Nneurons', default=10, type=int, help='Number of Neurons in each hidden layer')
  parser.add_argument('-M', '--Nlayers', default=2, type=int, help='Number of hidden layers')
  parser.add_argument('-NCPU', default=10, type=int, help='Number of Processors for TNC')
  parser.add_argument('--train_frac', default=.8, type=float, help='Fraction of LUT which is used to train network')

  args = parser.parse_args()

  print('Converting LUT in file   :::   {0:}  '.format(args.LUTfiles ))
  print('             coeffmode   :::   {0:}  '.format(args.coeffmode))
  print('             ANN file    :::   {0:}\n'.format(args.ANNname  ))
  print('             Nneurons    :::   {0:}\n'.format(args.Nneurons ))
  print('             Nlayers     :::   {0:}\n'.format(args.Nlayers  ))
  print('             NCPU        :::   {0:}\n'.format(args.NCPU     ))
  print('             Train fract :::   {0:}\n'.format(args.train_frac))

  try:
      train_dataset = [[],[],[],[]]
      for f in args.LUTfiles:
          print('Loading LUT file: ',f)
          tdata, pspace = create_training_dataset(args.coeffmode, f, training_fraction=args.train_frac)
          for i,v in enumerate(tdata):
            train_dataset[i].append(v)

      for i,v in enumerate(train_dataset):
          train_dataset[i] = np.vstack(v)

  except Exception as e:
      print('Error occured when creating training Datasets :',e)
      sys.exit(-1)

  try:
      train_ANN(args.ANNname, args.Nlayers, args.Nneurons, train_dataset, pspace, ncpu=args.NCPU)
  except Exception as e:
      print('Error occured when training network :',e)
      sys.exit(-1)
