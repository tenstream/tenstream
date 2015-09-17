#!/usr/bin/env python

import netCDF4 as NC
import numpy as np
import sys
import os
import shutil


def filename_from_url(url):
    import urlparse
    return os.path.basename(urlparse.urlsplit(url)[2])

def list_ftpdir(url):
    import ftplib
    import urlparse

    parsed_uri = urlparse.urlparse(url)
    domain = '{uri.scheme}://{uri.netloc}/'.format(uri=parsed_uri)    

    ftp = ftplib.FTP()
    ftp.connect( parsed_uri.netloc )
    ftp.login()

    ftp.cwd( os.path.dirname( parsed_uri.path ) )

    files = []

    try:
        files = ftp.nlst()
    except ftplib.error_perm, resp:
        if str(resp) == "550 No files found":
            print "No files in this directory"
        else:
            raise
    
#    print 'Server Repository has the following files in it: \n'
#    for f in files:
#        print f
#
#    print '\n Maybe you want one of those? \n'
    return files

def get_ftp_file(url):
    import ftplib
    import urlparse
    import tempfile

    parsed_uri = urlparse.urlparse(url)
    domain = '{uri.scheme}://{uri.netloc}/'.format(uri=parsed_uri)    

    ftp = ftplib.FTP()
    ftp.connect( parsed_uri.netloc )
    ftp.login()

    ftp.cwd( os.path.dirname( parsed_uri.path ) )

    fdst = tempfile.NamedTemporaryFile()
    print 'Downloading file from',url,' => ',fdst.name
    ftp.retrbinary('RETR %s' % os.path.basename( parsed_uri.path ), fdst.write)
    fdst.file.flush()

    ftp.quit()
    return fdst

def copy_nc_var(invar, Dout):
#    #Copy dimensions
#    for dname, the_dim in dsin.dimensions.iteritems():
#        print dname, len(the_dim)
#        dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    
    # Copy variables
    outVar = Dout.createVariable(invar.name, invar.datatype, invar.dimensions)

    # Copy variable attributes
    outVar.setncatts({k: invar.getncattr(k) for k in invar.ncattrs()})

    outVar[:] = invar[:]

def merge_nc_var(server,server_tol, local,local_tol):
    condition = server_tol[:] >= local_tol[:] # 1 if local better than server

    Nlocal  = np.sum(condition)
    Nserver = np.size(condition) - Nlocal
    print '    Using {0:} values from server and {1:} from local table'.format(Nserver,Nlocal)
    if Nlocal>0:
        print '      ATTENTION ::: If you have a lot of better coeffs than are available at the Server please send them to fabian@jakub.com'
    
    new     = np.choose( condition, [ server[:]    , local[:]      ] )
    new_tol = np.choose( condition, [ server_tol[:], local_tol[:]  ] )
    print 'max tol ( server,local ) ', np.max(server_tol), np.max(local_tol), ' :: new', np.max(new_tol)
    local[:] = new[:]
    local_tol[:] = new_tol[:]
    

def merge_LUT(LUT, serverLUT):
    print 'merging LUT {0:} ==> {1:}'.format(serverLUT,LUT)

    Dserver = NC.Dataset(serverLUT)
    Dlocal  = NC.Dataset(LUT, mode='a')

#    print 'Server LUTs',Dserver.variables.keys()
#    print 'Local  LUTs',Dlocal.variables.keys()
    print ''

    # Check if we can copy some tables
    for k in  Dserver.variables.keys():
        if k not in Dlocal.variables:
            print 'Found var which is not in local LUT... copying it over from server LUT :: ',k
            copy_nc_var(Dserver.variables[k], Dlocal)
        else:
            print 'already exists in local LUT :: ',k

    # Check if we should merge some values bc. server version has better tolerances
    for k in Dserver.variables.keys():
        if k.endswith('S') or k.endswith('T'):

            print '\nMerging NC Variable',k, k+'_tol'
            server     = Dserver.variables[k]
            server_tol = Dserver.variables[k+'_tol']
            local      = Dlocal.variables [k]
            local_tol  = Dlocal.variables [k+'_tol']
            merge_nc_var(server,server_tol, local,local_tol)


    Dserver.close()
    Dlocal.close()
    
    
def update_LUT(LUTpath, LUTserver):

    path,LUTname = os.path.split(LUTpath)

    try: # download LUT from ftp
        url = 'ftp://'+LUTserver+'/'+LUTname
        serverLUTfile = get_ftp_file(url)

    except Exception,e:
        print 'Error occured when we tried download LUT with name:',LUTname
        print e

        availfiles = list_ftpdir(url)
        print 'Server Repository has the following files in it: \n'
        for f in availfiles:
            print f
    
        print '\n Maybe you want one of those? \n'
        return

        
    if os.path.isfile(LUTpath): # merge LUT's
        merge_LUT(LUTpath, serverLUTfile.name)

    else: # just copy serverLUT to destination
        print 'Just downloaded LUT({0:}) and save to {1:}'.format(serverLUTfile.name,LUTpath)
        shutil.copyfile(serverLUTfile.name, LUTpath)
    
    serverLUTfile.close() # finally close temp file and with that, delete it

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Input parameters for LUT web downloader')
    parser.add_argument('LUTpath', type=str, help='LUT file which should be loaded/merged')
    parser.add_argument('-LUTserver', type=str, help='LUTserver which should be queried for LUTs', default='ftp.meteo.physik.uni-muenchen.de/public/TenStream_LUT/')
    parser.add_argument('-list', action='store_const',const=True, help='list all available files at the server', default=False)
    args = parser.parse_args()

    if args.list:
      availfiles = list_ftpdir('ftp://'+args.LUTserver)
      for f in availfiles:
        print f
      sys.exit(0)


    update_LUT(args.LUTpath, args.LUTserver)

