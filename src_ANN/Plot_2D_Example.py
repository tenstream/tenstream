import netCDF4 as nc
import numpy as np
import ffnet as ff
import matplotlib.pyplot as plt

LUT_file = '/home/users/max.eckl/scratch/LUT/LUT_8_10.direct.dx67.dy67.pspace.dz20.kabs20.ksca20.g3.phi0.theta0.delta_T_1.000.nc'
#LUT_file = '/home/users/max.eckl/scratch/LUT/LUT_8_10.diffuse.dx67.dy67.pspace.dz20.kabs20.ksca20.g3.delta_T_1.000.nc'

#net_file = '/home/users/max.eckl/scratch/ANN/dx67dy67/index/dir2dir/dir2dir_31_.net'
net_file = '/home/users/max.eckl/scratch/ANN/dx67dy67/index/dir2diff/dir2diff_10_.net'
#net_file = '/usr/users/max.eckl/scratch/ANN/dx67dy67/index/diff2diff/bak/diff2diff_13_.net'

net_res  = 51

dim      = 'kabs'

const_val    = { 
                  'dz'   :   9   ,
                  'kabs' :   0   ,
                  'ksca' :   0   ,
                  'g'    :   0   , 
                  'phi'  :   0   ,
                  'theta':   0   
                }


dir2dir=True



nplotx = 3
nploty = 3








def Get_Dim ( LUT_file, *key ):
   
    dim = {}
    for xkey in key:
        ind = LUT_file.rfind( xkey ) + len(xkey); num=0
        while LUT_file[ind+num]!='.': num=num+1
        dim[xkey] = int(LUT_file[ind:ind+num])

    return dim


###################################################################################################################
#                                                                                                                 #
#                                           load LUT data from nc file                                            #
#                                                                                                                 #
###################################################################################################################

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
    LUT_names.extend(['T'  ,'T_tol'])
    var      .extend(['phi','theta'])
    LUT_prefix = LUT_prefix+'.phi{0:d}.theta{1:d}'.format(res['phi'],res['theta'])
    T, T_tol = [], []

data = [LUT_data.variables['{0:s}.dx{1:d}.dy{2:d}.pspace.{3:s}'.format(coeff_type,res['dx'],res['dy'],xvar)][:] for xvar in var]
LUT  = [LUT_data.variables['{0:s}.{1:s}.{2:s}'.format(coeff_type,LUT_prefix,LUT_name)][:] for LUT_name in LUT_names]

for  dz  in range( res[ 'dz' ] ):
    for kabs in range( res['kabs'] ):
        for ksca in range( res['ksca'] ):
            for  g   in range( res[ 'g'  ] ):

                S    .append ( LUT[0][g,ksca,kabs,dz] )
                S_tol.append ( LUT[1][g,ksca,kabs,dz] )
                if coeff_type=='diffuse':
                    index.append ( [dz,kabs,ksca,g] )
                    src  .append ( [data[0][dz], data[1][kabs], data[2][ksca], data[3][g]] )
                else:
                    index.append ( [dz,kabs,ksca,g,np.argwhere(float(res['phi'])==data[4]),np.argwhere(float(res['theta'])==data[5])] )
                    src  .append ( [data[0][dz], data[1][kabs], data[2][ksca], data[3][g], float(res['phi']), float(res['theta'])] )
                    T    .append ( LUT[2][g,ksca,kabs,dz] )
                    T_tol.append ( LUT[3][g,ksca,kabs,dz] )

LUT = {'index':np.array(index), 'S':np.array(S), 'S_tol':np.array(S_tol)}; src=np.array(src)
for ind, xvar in enumerate(var): LUT[xvar]=np.array(data[ind])
if coeff_type=='direct': 
    LUT.update( {'T':np.array(T), 'T_tol':np.array(T_tol)} )
    app = np.append ( np.ones((len(src),1),dtype=float)*res['phi'], np.ones((len(src),1),dtype=float)*res['theta'], axis=1 )
    src = np.append ( src, app, axis=1 )
LUT.update( {'src':src} )

# result :: LUT = {'src','ind','S','S_tol',[ 'T','T_tol' ]}  [] if "direct-LUT-file"

###############################################################################################################################



dim_LUT  = LUT[dim]
dim_net  = np.linspace( np.min(dim_LUT), np.max(dim_LUT), num=net_res )
#TODO: ind<->val conversion
ind_net  = np.linspace( 0, len(dim_LUT)-1, num=net_res )
val, ind = [], []
for x in var:
    if x in LUT.keys():
        if x==dim:
            val.append(dim_net)
            ind.append(ind_net)
        else:
            val.append( np.ones((len(dim_net),),dtype=float) * LUT[x][const_val[x]] )
            ind.append( np.ones((len(dim_net),),dtype=float) * float(const_val[x]))

val, ind = np.array(val).T, np.array(ind).T



net = ff.loadnet(net_file)
net_out = net(ind)
coef_name = 'T' if dir2dir else 'S'
trg = LUT[coef_name]       .reshape( (len(LUT['dz']),len(LUT['kabs']),len(LUT['ksca']),len(LUT['g']),LUT[coef_name].shape[1]) )
err = LUT[coef_name+'_tol'].reshape( (len(LUT['dz']),len(LUT['kabs']),len(LUT['ksca']),len(LUT['g']),LUT[coef_name].shape[1]) )

if dim=='dz':
    trg = trg[:,const_val['kabs'],const_val['ksca'],const_val['g'],:]
    err = err[:,const_val['kabs'],const_val['ksca'],const_val['g'],:]
elif dim=='kabs':
    trg = trg[const_val['dz'],:,const_val['ksca'],const_val['g'],:]
    err = err[const_val['dz'],:,const_val['ksca'],const_val['g'],:]
elif dim=='ksca':
    trg = trg[const_val['dz'],const_val['kabs'],:,const_val['g'],:]
    err = err[const_val['dz'],const_val['kabs'],:,const_val['g'],:]
elif dim=='g':
    trg = trg[const_val['dz'],const_val['kabs'],const_val['ksca'],:,:]
    err = err[const_val['dz'],const_val['kabs'],const_val['ksca'],:,:]
else:
    raise ValueError('try another dim')

dim_ind = var.index(dim)

fig, ax = plt.subplots(nplotx,nploty,sharex=True,squeeze=False)
for i,xax in enumerate(ax):
    for j,yax in enumerate(xax):
        yax.plot(ind[:,dim_ind],net_out[:,i*nplotx+j],color='red',marker='D',label='ANN')
        yax.errorbar(np.arange(len(trg[:,i*nplotx+j]),dtype=float),trg[:,i*nplotx+j],yerr=err[:,i*nplotx+j]*10.0,color='blue',marker='s',label='LUT')
        yax.legend(loc='best'); yax.grid()
        yax.set_xlabel(dim,visible=True if i==2 else False); yax.set_ylabel('coeff {}'.format(i*3+j))
        yax.set_xlim(np.min(ind[:,dim_ind]),np.max(ind[:,dim_ind]))
        yax.set_xticks([0,4,8,12,16,19])
        yax.set_xticklabels(['{0:9.2e}'.format(LUT[dim][0]),'{0:9.2e}'.format(LUT[dim][4]),'{0:9.2e}'.format(LUT[dim][8]),'{0:9.2e}'.format(LUT[dim][12]),'{0:9.2e}'.format(LUT[dim][16]),'{0:9.2e}'.format(LUT[dim][19])],rotation=45)
       
title = [val[const_val[x],i] if i!=dim_ind else 'var' for i,x in enumerate(var)]
if coeff_type=='direct':
    fig.suptitle('dz :: {};  kabs :: {};  ksca :: {};  g :: {};  phi :: {};  theta :: {}'.format(*title))
else:
    fig.suptitle('dz :: {};  kabs :: {};  ksca :: {};  g :: {}'.format(*title))
fig.show()
