import numpy as np

def params() :
    '''
    Define the order of the parameter arrays, with the associated FERRE names, tag names, and flag names
    '''
    tagnames=np.array(['TEFF','LOGG','LOGVMICRO','M_H','C_M','N_M','ALPHA_M','LGVSINI','PARAM_O'])
    flagnames=np.array(['TEFF','LOGG','VMICRO','M_H','C_M','N_M','ALPHA_M','VSINI','O'])
    params=np.array(['TEFF','LOGG','LOG10VDOP','METALS','C','N','O Mg Si S Ca Ti','LGVSINI','O'])
    return params,tagnames,flagnames

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

       

