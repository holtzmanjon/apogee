"""
Routines for handling ASPCAP bitmasks
"""

import numpy as np

def getaspcapflags() :
    """
    get the ASPCAP flag names, descriptions, and warn/bad classification
    """

    flag=(['TEFF_WARN','LOGG_WARN','VMICRO_WARN','M_H_WARN','ALPHA_M_WARN','C_M_WARN','N_M_WARN','STAR_WARN',
          'CHI2_WARN','COLORTE_WARN','ROTATION_WARN','SN_WARN','SPEC_HOLE_WARN','ATMOS_HOLE_WARN','VSINI_WARN','',
          'TEFF_BAD','LOGG_BAD','VMICRO_BAD','M_H_BAD','ALPHA_M_BAD','C_M_BAD','N_M_BAD','STAR_BAD',
          'CHI2_BAD','COLORTE_BAD','ROTATION_BAD','SN_BAD','SPEC_HOLE_BAD','ATMOS_HOLE_BAD','VSINI_BAD','NO_ASPCAP_RESULT'])
    badflag=([2,2,0,0,0,0,0,2,
              2,2,2,2,2,2,0,0,
              1,1,0,0,0,0,0,1,
              1,1,1,1,1,2,0,1])

    descrip=([
     'WARNING on effective temperature (see PARAMFLAG[0] for details) ',
     'WARNING on log g (see PARAMFLAG[1] for details) ',
     'WARNING on vmicro (see PARAMFLAG[2] for details) ',
     'WARNING on metals (see PARAMFLAG[3] for details) ',
     'WARNING on [alpha/M] (see PARAMFLAG[4] for details) ',
     'WARNING on [C/M] (see PARAMFLAG[5] for details) ',
     'WARNING on [N/M] (see PARAMFLAG[6] for details) ',
     'WARNING overall for star: set if any of TEFF, LOGG, CHI2, COLORTE, ROTATION, SN warn are set ',
     'high chi^2 (> 2*median at ASPCAP temperature (WARN)',
     'effective temperature more than 500K from photometric temperature for dereddened color (WARN)',
     'Spectrum has broad lines, with possible bad effects: FWHM of cross-correlation of spectrum with best RV template relative to auto-correltion of template > 1.5 (WARN)',
     'S/N<70 (WARN)',
     'Grid point within 2 grid steps of hole-filled synthesis ',
     'Grid point within 2 grid steps of hole-filled atmosphere ',
     ' ',
     ' ',
     'BAD effective temperature (see PARAMFLAG[0] for details) ',
     'BAD log g (see PARAMFLAG[1] for details) ',
     'BAD vmicro (see PARAMFLAG[2] for details) ',
     'BAD metals (see PARAMFLAG[3] for details) ',
     'BAD [alpha/M] (see PARAMFLAG[4] for details) ',
     'BAD [C/M] (see PARAMFLAG[5] for details) ',
     'BAD [N/M] (see PARAMFLAG[6] for details) ',
     'BAD overall for star: set if any of TEFF, LOGG, CHI2, COLORTE, ROTATION, SN error are set, or any GRIDEDGE_BAD ',
     'high chi^2 (> 5*median at ASPCAP temperature (BAD)',
     'effective temperature more than 1000K from photometric temperature for dereddened color (BAD)',
     'Spectrum has broad lines, with possible bad effects: FWHM of cross-correlation of spectrum with best RV template relative to auto-correltion of template > 2 (BAD)',
     'S/N<50 (BAD)',
     'Grid point within 1 grid steps of hole-filled synthesis ',
     'Grid point within 1 grid steps of hole-filled atmosphere ',
     ' ',
     ' '
     ])
    names=[('flag','S24'),('badflag','i4'),('descrip','S80')]
    tmp=np.zeros(len(flag),dtype=names)
    tmp['flag']=flag 
    tmp['badflag']=badflag 
    tmp['descrip']=descrip 
    return tmp

def getstarflags() :
    """
    get the STARFLAG names, descriptions, and warn/bad classification
    """

    flag=(['BAD_PIXELS','COMMISSIONING','BRIGHT_NEIGHBOR','VERY_BRIGHT_NEIGHBOR','LOW_SNR','','','',
          '','PERSIST_HIGH','PERSIST_MED','PERSIST_LOW','PERSIST_JUMP_POS','PERSIST_JUMP_NEG','','',
          'SUSPECT_RV_COMBINATION','SUSPECT_BROAD_LINES','BAD_RV_COMBINATION','','','','','',
          '','','','','','','',''])
    badflag=([1,0,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,1,0,0,0,0,0,
             0,0,0,0,0,0,0,0]   )

    descrip=([
     'Spectrum has many bad pixels (>20%):  BAD',                                                         
     'Commissioning data (MJD<55761), non-standard configuration, poor LSF: WARN',                       
     'Star has neighbor more than 10 times brighter: WARN',
     'Star has neighbor more than 100 times brighter: BAD',
     'Spectrum has low S/N (S/N<5)',                                                                    
     '',
     '',
     '',
     '',
     'Spectrum has significant number (>20%) of pixels in high persistence region: WARN',               
     'Spectrum has significant number (>20%) of pixels in medium persistence region: WARN',
     'Spectrum has significant number (>20%) of pixels in low persistence region: WARN',
     'Spectrum show obvious positive jump in blue chip: WARN',
     'Spectrum show obvious negative jump in blue chip: WARN',                                         
     '',
     '',
     'RVs from synthetic template differ significantly (~2 km/s) from those from combined template: WARN', 
     'Cross-correlation peak with template significantly broader than autocorrelation of template: WARN',
     'RVs from synthetic template differ very significatly (~10 km/s) from those from combined template: BAD',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     ''
     ])
    names=[('flag','S24'),('badflag','i4'),('descrip','S80')]
    tmp=np.zeros(len(flag),dtype=names)
    tmp['flag']=flag 
    tmp['badflag']=badflag 
    tmp['descrip']=descrip 
    return tmp

def getpixmask() :

    flag=(['BADPIX','CRPIX','SATPIX','UNFIXABLE','BADDARK','BADFLAT','BADERR','NOSKY',
          'LITTROW_GHOST','PERSIST_HIGH','PERSIST_MED','PERSIST_LOW','SIG_SKYLINE','SIG_TELLURIC','NOT_ENOUGH_PSF',''])

    badflag=([1,1,1,1,1,1,1,1,
             0,0,0,0,0,0,1,0])

    maskcontrib=([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
                 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0])

    descrip=([
     'Pixel marked as BAD in bad pixel mask or from strong persistence jump',
     'Pixel marked as cosmic ray in ap3d',
     'Pixel marked as saturated in ap3d',
     'Pixel marked as unfixable in ap3d',
     'Pixel marked as bad as determined from dark frame',
     'Pixel marked as bad as determined from flat frame',
     'Pixel set to have very high error (not used)',
     'No sky available for this pixel from sky fibers',
     'Pixel falls in Littrow ghost, may be affected',
     'Pixel falls in high persistence region, may be affected',
     'Pixel falls in medium persistence region, may be affected',
     'Pixel falls in low persistence region, may be affected',
     'Pixel falls near sky line that has significant flux compared with object',
     'Pixel falls near telluric line that has significant absorption',
     'Less than 50 percent PSF in good pixels',
     ''
    ])
    names=[('flag','S24'),('badflag','i4'),('descrip','S80'),('maskcontrib','f4')]
    tmp=np.zeros(len(flag),dtype=names)
    tmp['flag']=flag 
    tmp['badflag']=badflag 
    tmp['descrip']=descrip 
    tmp['maskcontrib']=maskcontrib 
    return tmp


# old name compatibilities
def aspcapflagval(flag) :
    return val('ASPCAPFLAG',flag)
def badaspcapflag() :
    return bad('ASPCAPFLAG')
def warnaspcapflag() :
    return warn('ASPCAPFLAG')
def aspcapflag(mask,type=0) :
    return flag('ASPCAPFLAG',mask,type=0)
def persist() :
    return val('STARFLAG','PERSIST_HIGH') | val('STARFLAG','PERSIST_MED') | val('STARFLAG','PERSIST_LOW') | val('STARFLAG','PERSIST_JUMP_POS') | val('STARFLAG','PERSIST_JUMP_NEG') 

def getflags(bitmask) :
    if bitmask.upper() == 'ASPCAPFLAG' :
        flags = getaspcapflags()
    elif bitmask.upper() == 'STARFLAG' :
        flags = getstarflags()
    elif bitmask.upper() == 'PIXMASK' :
        flags = getpixmask()
    return flags

def val(bitmask,flag) :
    """
    Get the numerical bit value of a given character ASPCAP flag
    """
    flags=getflags(bitmask)
    j=np.where(flags['flag'] == flag.strip())[0]
    if len(j) > 0 :
        bitval=2**j[0] 
    else :
        bitval=0
        print('WARNING: undefined mask: ',flag)
    return bitval


def bad(bitmask) :
    """
    Return bitmask of values that indicate BAD in input bitmask
    """
    flags=getflags(bitmask)
    bad=0
    i=0
    for flag in flags['badflag'] :
        if flag == 1 :
            bad=bad | 2**i
        i+=1
    return bad


def warn(bitmask) :
    """
    Return bitmask of values that indicate WARN or BAD in input bitmask
    """
    flags=getflags(bitmask)
    bad=0
    i=0
    for flag in flags['badflag'] :
        if flag >= 1 :
            bad=bad | 2**i
        i+=1
    return bad

def flag(bitmask,mask,type=0) :
    """
    Return names for all bits that are set in input mask for input bitmask
    """
    flags=getflags(bitmask)
    strflag=''
    ibit = 0
    for name in flags['flag'] :
        if ( mask & 2**ibit ) > 0 and ( type == 0 or flag['badflag'] == type ) :
          strflag = strflag + name +','
        ibit+=1

    return strflag.strip(',')

