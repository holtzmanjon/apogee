import numpy as np
import matplotlib.pyplot as plt
import os
import pdb
import scipy.ndimage.filters
from astropy.io import fits
from astropy.io import ascii
from holtz.tools import struct
from holtz.tools import plots
from holtz.apogee import apload

def params() :
    '''
    Define the order of the parameter arrays, with the associated FERRE names, tag names, and flag names
    '''
    tagnames=np.array(['TEFF','LOGG','LOGVMICRO','M_H','C_M','N_M','ALPHA_M','LGVSINI','PARAM_O'])
    flagnames=np.array(['TEFF','LOGG','VMICRO','M_H','C_M','N_M','ALPHA_M','VSINI','O'])
    params=np.array(['TEFF','LOGG','LOG10VDOP','METALS','C','N','O Mg Si S Ca Ti','LGVSINI','O'])
    return params,tagnames,flagnames

def elems(nelem=0) :
    '''
    Define the order of the element arrays
    '''

    elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
    #return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni','Nd']
    #return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni']
    elemtoh=[0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] #0,0,0,0,0,0,0,0,0]
    tagnames=[]
    elemfitnames=[]
    for i in range(len(elems) ) :
        if elemtoh[i] :
            tagnames.append(elems[i]+'_Fe')
            elemfitnames.append('['+elems[i]+'/Fe]')
        else :
            tagnames.append(elems[i]+'_M')
            elemfitnames.append('['+elems[i]+'/M]')
    return elems,elemtoh,tagnames,elemfitnames

def readstars(starlist,libpar) :
    '''
    Runs stars in starlist through FERRE using libpar
    '''
    for star in starlist :
        spec,err=readstar(star)
        cont=cont_normalize(spec)

def ferre(spec,libpar,te0=None,logg0=None,mh0=None) :
    '''
    Runs FERRE for a set of input spectra
    '''
    np.savetxt(out+'.frd',spec)
    np.savetxt(out+'.err',spec)
    np.savetxt(out+'.con',spec)

    # write ferre control file
    ferre.writeferre()

    out=ferre.readferre(name)

def elemsens(els=None,plot=None,ylim=[0.1,-0.3],teff=4750,logg=2.,feh=-1.,smooth=None) :
    '''
    Returns and optionally plots wavelength sensitivity to changes in elemental abundances for specified elements from MOOG mini-elem grid
    '''
    elem=fits.open(os.environ['APOGEE_REDUX']+'/speclib/moog/elemsens.fits')
    if els is None :
        els = elems()[0]
    elif type(els) == str :
        els = [els]
    wave=[]
    out=[]
    for el in els :
        for i in range(1,25) :
            card='HDU{:02d}'.format(i)
            try :
              if elem[0].header[card].strip().upper() == el.strip().upper() :
                it=int(round((teff-elem[i].header['CRVAL2'])/elem[i].header['CDELT2']))
                ig=int(round((logg-elem[i].header['CRVAL3'])/elem[i].header['CDELT3']))
                ife=int(round((feh-elem[i].header['CRVAL4'])/elem[i].header['CDELT4']))
                diff=elem[i].data[ife,ig,it,:]
                if smooth is not None:
                    diff=scipy.ndimage.filters.gaussian_filter(diff,smooth)
                wave=elem[i].header['CRVAL1']+np.arange(elem[i].header['NAXIS1'])*elem[i].header['CDELT1']
                if plot is not None:
                    #plot.plot(wave,diff,color='g')
                    plot.plot(wave,diff)
                    plot.set_ylim(ylim[0],ylim[1])
                    plt.draw()
                out.append(diff)
            except: pass
    if len(out) == 1 :
        return wave, out[0]
    else :
        return wave, out

def sensplot(ax=None,offset=0) :
    if ax is None :
        fig,ax=plots.multi(1,2,hspace=0.001,sharex=True)
    els=['O','Mg','Si','S','Ca','Ti','Na','Al','K','P']
    cols=['r','g','b','c','y','m','r','g','b','c']
    ls=['-','-','-','-','-','-',':',':',':',':']
    for i in range(len(els)) :
        w,s=elemsens(els=els[i])
        plots.plotl(ax[0],w,s+offset,label=els[i],color=cols[i],ls=ls[i])
    ax[0].legend(fontsize='small')

    elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
    els=['V','Cr','Mn','Co','Ni','Cu','Ge','Ce','Rb','Nd']
    cols=['r','g','b','c','y','m','r','g','b','c']
    ls=['-','-','-','-','-','-',':',':',':',':']
    for i in range(len(els)) :
        w,s=elemsens(els=els[i])
        plots.plotl(ax[1],w,s+offset,label=els[i],color=cols[i],ls=ls[i])
    ax[1].legend(fontsize='small')
    #elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
    

def data(str,loc=None) :       
    '''
    Add apogeeObject data to structure
    '''
    add=np.empty(1,dtype=[('RA','f4'),('DEC','f4'),('J','f4'),('H','f4'),('K','f4'),('AK_TARG','f4'),('SFD_EBV','f4')])
    new=struct.add_cols(str,add)
    for i in range(len(str)) :
        name = str['APOGEE_ID'][i]
        print i, name
        apogee_id = name.split('.')[0].split('-')[2]
        if loc is None :
            loc = name.split('.')[1].split('_')[0]
        s=apload.apStar(loc,apogee_id)
        field=s[0].header['FIELD']
        try :
            obj=fits.open(os.environ['APOGEE_TARGET']+'/apogeeObject/apogeeObject_'+field+'.fits')[1].data
        except :
            obj=fits.open(os.environ['APOGEE_TARGET']+'/apogee2object/apogee2object_'+field+'.fits')[1].data
        j=np.where(obj['APOGEE_ID'] == apogee_id)[0]
        for card in add.dtype.names :
            new[card][i]=obj[card][j]
    return new

def periodic(n) :
    '''
    Routine to get element name / atomic number conversion
    '''
    elem=np.array(['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd'])
    if type(n) is str :
        j=np.where(elem == n)[0]
        return j+1
    else :
        if n == 0 :
            return ''    
        else :
            return elem[n-1]

def solar(el) :
    sunabund_2007  =  np.array([
       12.00, 10.93,  1.05,  1.38,  2.70,  8.39,  7.78,  8.66,  4.56,   #  1 -  9
        7.84,  6.17,  7.53,  6.37,  7.51,  5.36,  7.14,  5.50,  6.18,   # 10 - 18
        5.08,  6.31,  3.17,  4.90,  4.00,  5.64,  5.39,  7.45,  4.92,   # 19 - 27
        6.23,  4.21,  4.60,  2.88,  3.58,  2.29,  3.33,  2.56,  3.25,   # 28 - 36
        2.60,  2.92,  2.21,  2.58,  1.42,  1.92, -99.0,  1.84,  1.12,   # 37 - 45
        1.66,  0.94,  1.77,  1.60,  2.00,  1.00,  2.19,  1.51,  2.24,   # 46 - 54
        1.07,  2.17,  1.13,  1.70,  0.58,  1.45, -99.0,  1.00,  0.52,   # 55 - 63
        1.11,  0.28,  1.14,  0.51,  0.93,  0.00,  1.08,  0.06,  0.88,   # 64 - 72
       -0.17,  1.11,  0.23,  1.25,  1.38,  1.64,  1.01,  1.13,  0.90,   # 73 - 81
        2.00,  0.65, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.06,   # 82 - 90
       -99.0, -0.52 ])                                                  # 91 - 9
    n = periodic(el)
    return(sunabund_2007[n-1])

def rydberg(n1,n2) :
    '''
    Rydberg formula to give H wavelengths for transitions between levels n1 and n2
    '''
    r=1.0973731568e7
    me=9.109382e-31
    mprot=1.672621e-27
    rm=r/(1+me/mprot)

    l=rm*(1./n1**2-1./n2**2)
    w=1./l*1.e10
    return w

def hlines(plot=None,yloc=0.,n1=4,n2=range(11,22)) :
    '''
    Return approximate (Rydberg) location of H lines, defaulting to lines in APOGEE spectra
    '''
    h=[]
    for n in n2 :
        h.append(rydberg(n1,n))
    h=np.array(h)
    if plot is not None :
        plots.plotp(plot,h,h*0.+yloc)
    return h

def elemmask(el,maskdir='filters_26112015',plot=None,yr=[0,1]) :
    '''
    '''
    mask=np.loadtxt(os.getenv('SPECLIB_DIR')+'/lib/'+maskdir+'/'+el+'.filt') 
    wave=np.loadtxt(os.getenv('SPECLIB_DIR')+'/lib/'+maskdir+'/wave.dat')
    if plot is not None :
        plots.plotl(plot,wave,mask,yr=yr)
    return wave,mask
