from sdss.apogee import apload
from sdss.apogee import apselect
from astropy.io import fits
import pdb
import numpy as np
import os

def main() :
    '''
    
    '''

    apload.dr13()
    a=fits.open('APOGEE_GES_overlap_calib_DR13.fits')[1].data

    for i in range(len(a)) :
        if a['LOCATION_ID'][i] == 1 :
          apload.aspcapStar(a['FIELD'][i],a['REDUCTION_ID'][i])
          apload.apStar1m(a['FIELD'][i],a['REDUCTION_ID'][i])
        else :
          apload.aspcapStar(a['LOCATION_ID'][i],a['REDUCTION_ID'][i])
          apload.apStar(a['LOCATION_ID'][i],a['REDUCTION_ID'][i])
