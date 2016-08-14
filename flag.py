"""
Routines for handling ASPCAP flags
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


def aspcapflagval(flag) :
    """
    Get the numerical bit value of a given character ASPCAP flag
    """
    flags = getaspcapflags()
    j=np.where(flags['flag'] == flag.strip())[0]
    if len(j) > 0 :
        val=2**j[0] 
    else :
        val=0
        print('WARNING: undefined mask: ',flag)
    return val


def badaspcapflag() :
    """
    Return bitmask of values that indicate BAD ASPCAP
    """
    flags = getaspcapflags()
    bad=0
    i=0
    for flag in flags['badflag'] :
        if flag == 1 :
            bad=bad | 2**i
        i+=1
    return bad


def warnaspcapflag() :
    """
    Return bitmask of values that indicate WARN or BAD ASPCAP
    """
    flags = getaspcapflags()
    bad=0
    i=0
    for flag in flags['badflag'] :
        if flag >= 1 :
            bad=bad | 2**i
        i+=1
    return bad

def aspcapflag(mask,type=0) :

    strflag=''
    flag = getaspcapflags()
    ibit = 0
    for name in flag['flag'] :
        if ( mask & 2**ibit ) > 0 and ( type == 0 or flag['badflag'] == type ) :
          strflag = strflag + name +','
        ibit+=1

    return strflag.strip(',')

