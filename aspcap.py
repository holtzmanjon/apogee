def params() :
    '''
    Define the order of the parameter arrays, with the associated FERRE names, tag names, and flag names
    '''
    tagnames=['TEFF','LOGG','LOGVMICRO','M_H','C_M','N_M','ALPHA_M','LGVSINI','PARAM_O']
    flagnames=['TEFF','LOGG','VMICRO','M_H','C_M','N_M','ALPHA_M','VSINI','O']
    params=['TEFF','LOGG','LOG10VDOP','METALS','C','N','O Mg Si S Ca Ti','LGVSINI','O']
    return params,tagnames,flagnames

