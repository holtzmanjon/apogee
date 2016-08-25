###############################################################################
#
#   apogee.download: download APOGEE data files
#
###############################################################################
import os
import sys
import shutil
import tempfile
import subprocess
import numpy
import pdb

from sdss_access.path import path
from sdss_access.sync.http import HttpAccess

_MAX_NTRIES= 2
_ERASESTR= "                                                                                "

def allfile(root,dr=None,apred=None,apstar=None,aspcap=None,results=None,location=None,obj=None,plate=None,mjd=None,num=None,telescope='apo25m',fiber=None,chips=False):
    """ download the allStar file """
    sdss_path=path.Path()
    http_access=HttpAccess(verbose=True)
    http_access.remote()

    if chips == False :
        # First make sure the file doesn't exist locally
        filePath = sdss_path.full(root,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber)
        downloadPath = sdss_path.url(root,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber)
        if os.path.exists(filePath) is False: 
            http_access.get(root,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber)
        return filePath
    else :
        for chip in ['a','b','c'] :
            filePath = sdss_path.full(root,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber,
                chip=chip)
            if os.path.exists(filePath) is False: 
                http_access.get(root,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                    location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber,
                chip=chip)
        return filePath.replace('-c','')

def _download_file(downloadPath,filePath,dr,verbose=True,spider=False):

    print('downloading: ', downloadPath, ' to: ', filePath)
    #return

    try:
        # make all intermediate directories
        os.makedirs(os.path.dirname(filePath)) 
    except OSError: pass
    # Safe way of downloading
    downloading= True
    interrupted= False
    file, tmp_savefilename= tempfile.mkstemp()
    os.close(file) #Easier this way
    ntries= 0
    while downloading:
        try:
            cmd= ['wget','%s' % downloadPath,
                  '-O','%s' % tmp_savefilename]
            if not verbose: cmd.append('-q')
            if spider: cmd.append('--spider')
            if verbose: print('cmd: ',cmd)
            subprocess.check_call(cmd)
            if not spider: shutil.move(tmp_savefilename,filePath)
            downloading= False
            if interrupted:
                raise KeyboardInterrupt
        except subprocess.CalledProcessError:
            if not downloading: #Assume KeyboardInterrupt
                raise
            elif ntries > _MAX_NTRIES:
                raise IOError('File %s does not appear to exist on the server ...' % (os.path.basename(filePath)))
            sys.stdout.write('\r'+"KeyboardInterrupt ignored while downloading ...\r")
            sys.stdout.flush()
            os.remove(tmp_savefilename)
            interrupted= True
            ntries+= 1
        finally:
            if os.path.exists(tmp_savefilename):
                os.remove(tmp_savefilename)   
            os.chmod(filePath,0664)
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()        
    return None
