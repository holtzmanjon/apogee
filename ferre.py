import numpy as np
import pdb
from holtz.tools import plots
from holtz.tools import match
from holtz.apogee import aspcap

def writespec(name,data) :
    '''
    Writes FERRE "spectrum" file with input data, one line per star
    '''
    f=open(name,'w')
    for spec in data :
        for pix in np.arange(spec.shape[0]) :
            f.write('{:12.6f}'.format(spec[pix]))
        f.write('\n')
    f.close()
    return 


def writeferre() :
    '''
    Writes FERRE control file
    '''

def writeipf(name,libfile,stars,param=None) :
    '''
    Writes FERRE input file
    '''
    # get library headers and load wavelength array
    libhead0, libhead=rdlibhead(libfile)

    params=aspcap.params()[0]
    nparams=libhead0['N_OF_DIM']
    # get the index numbers in input parameter array for correct library parameter order
    index=np.zeros(nparams,dtype=int)
    for i in range(nparams) :
        index[i] = np.where(params == libhead0['LABEL'][i])[0]
    # if input parameters aren't specified, use zeros
    if param is None :
        param=np.zeros(len(params))
    # write the IPF file
    f=open(name+'.ipf','w')
    i=0
    for star in stars :
        f.write('{:<40s}'.format(star))
        for ipar in range(nparams) : 
            f.write('{:12.3f}'.format(param[i][index[ipar]]))
        f.write('\n')
        i+=1
    f.close()
        

def read(name,libfile) :
    '''
    Read all of the FERRE files associated with a FERRE run
    '''
    # get library headers and load wavelength array
    libhead0, libhead=rdlibhead(libfile)
    wave=[]
    for ichip in range(len(libhead)) :
        wave.extend(libhead[ichip]['WAVE'][0] + np.arange(libhead[ichip]['NPIX'])*libhead[ichip]['WAVE'][1])
    wave=10.**np.array(wave)
    nwave=len(wave)
    nparam=libhead0['N_OF_DIM']

    # input ipf and output spm files, match object names
    ipfobj,ipf=readferredata(name+'.ipf')
    spmobj,spm=readferredata(name+'.spm')
    i1,i2=match.match(ipfobj,spmobj)
    nobj=len(ipfobj)
    param=spm[:,0:nparam]
    paramerr=spm[:,nparam:2*nparam]
    chi2=10.**spm[:,2*nparam+2]
    covar=spm[:,2*nparam+3:2*nparam+3+nparam**2]
    covar=np.reshape(covar,(nobj,nparam,nparam))

    # load param array
    params=aspcap.params()[0]
    ntotparams=len(params)
    index=np.zeros(ntotparams,dtype=int)
    a=np.zeros(nobj, dtype=[('APOGEE_ID','S100'),
                              ('FPARAM','f4',(ntotparams)),
                              ('FPARAM_COV','f4',(ntotparams,ntotparams)),
                              ('PARAM_CHI2','f4')])
    a['APOGEE_ID']=ipfobj
    for i in range(ntotparams) :
        try :
            index[i] = np.where(libhead0['LABEL'] == params[i])[0]
            a['FPARAM'][:,i] = param[:,index[i]]
            for j in range(ntotparams) :
                a['FPARAM_COV'][:,i,j]=covar[:,index[i],index[j]]
        except :
            index[i] = -1
    a['PARAM_CHI2']=chi2

    # put it all into a structured array
    form='{:d}f4'.format(nwave)
    sform='{:d}f4'.format(spm.shape[1])
    out=np.empty(nobj, dtype=[('obj','S24'),('spm',sform),('obs',form),('err',form),('mdl',form),('chi2',form)])
    out['obj']=ipfobj
    out['spm'][i1,:]=spm[i2,:]
    out['obs']=readspec(name+'.frd')[i2,:]
    out['err']=readspec(name+'.err')[i2,:]
    out['mdl']=readspec(name+'.mdl')[i2,:]
    out['chi2']=(out['obs']-out['mdl'])**2/out['err']**2

    return a,out,wave

def readspec(name) :
    '''
    Read a single file with FERRE-format spectra, and return as 2D array [nspec,nwave]
    '''
    f=open(name)
    data=[]
    for line in f :
        spec=np.array(line.split())
        spec=spec.astype(float)
        data.append(spec)
    return np.array(data)

def readferredata(name) :
    '''
    Read a single file with FERRE-format data, and return as 2D array [nspec,nwave]
    '''
    f=open(name)
    alldata=[]
    allobj=[]
    for line in f :
        a=line.split()
        obj=np.array(a[0])
        data=np.array(a[1:len(a)])
        data=data.astype(float)
        alldata.append(data)
        allobj.append(obj)
    return np.array(allobj),np.array(alldata)


def rdsinglehead(f) :
    '''
    Read a FERRE library header into a dictionary
    '''
    dict={}
    for line in f:  
      # read until we hit a '/' in first character
      if line[1] is not '/' : 
        words=line.split()
        nwords=len(words)
        card=words[0].upper()
        if card[-1]=='=' : 
            # adjust if = is not separated from card name
            card = card[0:-1]
            nwords+=1
        if nwords>2 :
            # we have a value, set it to characters to right of = sign, not including newling
            value=line[line.find('=')+1:len(line)-1]
            if value.find("'") >= 0 :
                # we have a string, strip it
                val=value.replace("'"," ").strip()
                if card.find('LABEL') >= 0 :
                    ilab= int(card[card.find('(')+1:card.find(')')])
                    label[ilab-1] = val
            else :
                # we have numbers, put into array of correct type
                vals=value.split()
                n=len(vals)
                try:
                    val=np.array(vals).astype(int)
                except :
                    val=np.array(vals).astype(float)
                # if we just have one variable, we don't want an array
                if n == 1: val=val[0]
            if card == 'N_OF_DIM' : label=np.zeros(val,dtype='S24')
            if card.find('LABEL') < 0 :
                # add card, value to dictionary if not a LABEL() card
                dict[card]=val
      else:
          break
    try:
        dict['LABEL'] = label
    except:
        pass
    return dict

def rdlibhead(name) :
    '''
    Read a full FERRE library header with multi-extensions

    Returns:
       libstr, libstr : first header, then list of extension headers; headers returned as dictionaries
    '''
    try:
        f=open(name,'r')
    except:
        print('cannot open file',name)
        return
    libstr0=rdsinglehead(f)
    try:
        multi=libstr0['MULTI']
        libstr=[]
        for imulti in range(multi) :
            libstr.append(rdsinglehead(f))
    except:
        pass

    return libstr0,libstr

def plotspec(w,spec,n=0) :
    fig,ax=plots.multi(1,2)
    plots.plotl(ax[0],w,spec['obs'][n,:],yr=[0,1.3]) 
    plots.plotl(ax[0],w,spec['err'][n,:]) 
    plots.plotl(ax[0],w,spec['mdl'][n,:]) 
    plots.plotl(ax[1],w,spec['chi2'][n,:])
