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

from sdss.files import path

_DR10_URL= 'http://data.sdss3.org/sas/dr10/apogee'
_DR12_URL= 'http://data.sdss3.org/sas/dr12/apogee'
_PROPRIETARY_URL= 'https://data.sdss.org/sas/apogeework/apogee'
_MAX_NTRIES= 2
_ERASESTR= "                                                                                "

def allfile(root,dr=None,apred=None,apstar=None,aspcap=None,results=None,location=None,obj=None,plate=None,mjd=None,num=None,telescope='apo25m',fiber=None,chips=False):
    """ download the allStar file """
    sdss_path=path()
   
    # First make sure the file doesn't exist locally
    filePath, relPath = sdss_path.full_relative(root,apred=apred,apstar=apstar,aspcap=aspcap,results=results,location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber)
    print(filePath, relPath)
    if chips == False :
        if os.path.exists(filePath): return filePath

        # Create the remote file path
        downloadPath = ( _base_url(dr=dr)+'/spectro/redux/'+ relPath)
        _download_file(downloadPath,filePath,dr,verbose=True)
        return filePath
    else :
        for chip in ['a','b','c'] :
            file = filePath.replace(root,root+'-'+chip)
            if os.path.exists(file): 
                pass
            else :
                # Create the remote file path
                downloadPath = (_base_url(dr=dr)+'/spectro/redux/'+
                                 relPath.replace(root,root+'-'+chip) )
                _download_file(downloadPath,file,dr,verbose=True)
        return filePath

def _base_url(dr,rc=False):
    if dr is None : return _PROPRIETARY_URL

    if dr.upper() == 'DR10': return _DR10_URL
    elif dr.upper() == 'DR12': return _DR12_URL
    else: return _PROPRIETARY_URL

def _dr_string(dr):
    if dr == 'bosswork': return 'bosswork'
    else: return 'dr%s' % dr
 

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
