# routines for comparing gravities with asteroseismic sample

from sdss.apogee import apload
from sdss.apogee import apselect
from astropy.io import fits
from holtz.gal import isochrones
from holtz.tools import match
from holtz.tools import plots
from holtz.apogee import flag
import pdb
import matplotlib.pyplot as plt
import numpy as np
import os
import astropy

def dr13() :
    '''
    
    '''

    apload.dr13()
    dr13=apload.allStar()[1].data
    irfm=fits.open(os.environ['IDLWRAP_DIR']+'/data/irfm_temp.fits')[1].data
    i1,i2=match.match(dr13['APOGEE_ID'],irfm['2MASS ID'])
    pdb.set_trace()

    fig,ax=plots.multi(2,3,hspace=0.5,wspace=0.001)
    # diff color-coded by gravity as f([M/H])
    plots.plotc(ax[0,0],dr13['FPARAM'][i1,3],dr13['FPARAM'][i1,0]-irfm['IRFM TEFF'][i2],dr13['FPARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-300,300],xt='[M/H]',yt='ASPCAP-seismic Teff')
    plots.plotc(ax[0,1],dr13['PARAM'][i1,3],dr13['PARAM'][i1,0]-irfm['IRFM TEFF'][i2],dr13['PARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-300,300],xt='[M/H]',colorbar=True,zt='Teff')

    # diff color-coded by [M/H] as f(log g)
    plots.plotc(ax[1,0],dr13['FPARAM'][i1,1],dr13['FPARAM'][i1,0]-irfm['IRFM TEFF'][i2],dr13['FPARAM'][i1,3],xr=[0,5],zr=[-2.5,0.5],yr=[-300,300],zt='[M/H]',yt='ASPCAP-seismic Teff',xt='log g')
    plots.plotc(ax[1,1],dr13['PARAM'][i1,1],dr13['PARAM'][i1,0]-irfm['IRFM TEFF'][i2],dr13['PARAM'][i1,3],xr=[0,5],zr=[-2.5,0.5],yr=[-300,300],zt='[M/H]',colorbar=True,xt='log g')

