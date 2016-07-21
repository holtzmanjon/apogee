# routines for comparing gravities with photometric` sample

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

def bindata(xdata,ydata,bins) :
    mean=bins*0.
    for i in range(len(bins)-1) :
      j=np.where((xdata>bins[i]) & (xdata<bins[i+1]))[0]
      mean[i]=ydata[j].mean() 
    return mean

def dr13() :
    '''
    
    '''

    apload.dr13()
    dr13=apload.allStar()[1].data
    irfm=fits.open(os.environ['IDLWRAP_DIR']+'/data/irfm_temp.fits')[1].data
    i1,i2=match.match(np.chararray.strip(dr13['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID']))
    
    dr13=dr13[i1]
    irfm=irfm[i2]
    saga=np.where(irfm['SOURCE'] == 'SAGA')
    cl=np.where(irfm['SOURCE'] == 'CL')
    th=np.where(irfm['SOURCE'] == 'TH')
    sfd=np.where(irfm['SOURCE'] == 'SFD')

    fig,ax=plots.multi(2,2,hspace=0.001,wspace=0.001)
    bins=np.arange(-2.5,0.75,0.25)
    # diff color-coded by gravity as f([M/H])
    plots.plotc(ax[0,0],dr13['FPARAM'][saga,3],dr13['FPARAM'][saga,0]-irfm['IRFM TEFF'][saga],dr13['FPARAM'][saga,0],zr=[3500,5500],xr=[-2.5,0.5],yr=[-300,300],xt='[M/H]',yt='ASPCAP-photometric Teff')
    mean=bindata(dr13['FPARAM'][saga,3][0],(dr13['FPARAM'][saga,0]-irfm['IRFM TEFF'][saga])[0],bins)
    plots.plotp(ax[0,0],bins,mean,marker='o',size=40)
    ax[0,0].text(0.1,0.9,'SAGA',transform=ax[0,0].transAxes)

    plots.plotc(ax[0,1],dr13['FPARAM'][cl,3],dr13['FPARAM'][cl,0]-irfm['IRFM TEFF'][cl],dr13['FPARAM'][cl,0],zr=[3500,5500],xr=[-2.5,0.5],yr=[-300,300],xt='[M/H]')
    mean=bindata(dr13['FPARAM'][cl,3][0],(dr13['FPARAM'][cl,0]-irfm['IRFM TEFF'][cl])[0],bins)
    plots.plotp(ax[0,1],bins,mean,marker='o',size=40)
    ax[0,1].text(0.1,0.9,'CL',transform=ax[0,1].transAxes)

    plots.plotc(ax[1,0],dr13['FPARAM'][th,3],dr13['FPARAM'][th,0]-irfm['IRFM TEFF'][th],dr13['FPARAM'][th,0],zr=[3500,5500],xr=[-2.5,0.5],yr=[-300,300],xt='[M/H]',yt='ASPCAP-photometric Teff')
    mean=bindata(dr13['FPARAM'][th,3][0],(dr13['FPARAM'][th,0]-irfm['IRFM TEFF'][th])[0],bins)
    plots.plotp(ax[1,0],bins,mean,marker='o',size=40)
    ax[1,0].text(0.1,0.9,'TH',transform=ax[1,0].transAxes)

    plots.plotc(ax[1,1],dr13['FPARAM'][sfd,3],dr13['FPARAM'][sfd,0]-irfm['IRFM TEFF'][sfd],dr13['FPARAM'][sfd,0],zr=[3500,5500],xr=[-2.5,0.5],yr=[-300,300],xt='[M/H]')
    mean=bindata(dr13['FPARAM'][sfd,3][0],(dr13['FPARAM'][sfd,0]-irfm['IRFM TEFF'][sfd])[0],bins)
    plots.plotp(ax[1,1],bins,mean,marker='o',size=40)
    ax[1,1].text(0.1,0.9,'SFD',transform=ax[1,1].transAxes)

    fig.savefig('dteff_mh.jpg')
    pdb.set_trace()

    fig,ax=plots.multi(2,2,hspace=0.001,wspace=0.001)
    bins=np.arange(3500,5500,250)
    # diff color-coded by gravity as f([M/H])
    plots.plotc(ax[0,0],dr13['FPARAM'][saga,0],dr13['FPARAM'][saga,0]-irfm['IRFM TEFF'][saga],dr13['FPARAM'][saga,3],zr=[-2.5,0.5],xr=[6000,3000],yr=[-300,300],xt='Teff',yt='ASPCAP-photometric Teff')
    mean=bindata(dr13['FPARAM'][saga,0][0],(dr13['FPARAM'][saga,0]-irfm['IRFM TEFF'][saga])[0],bins)
    plots.plotp(ax[0,0],bins,mean,marker='o',size=40)
    ax[0,0].text(0.1,0.9,'SAGA',transform=ax[0,0].transAxes)

    plots.plotc(ax[0,1],dr13['FPARAM'][cl,0],dr13['FPARAM'][cl,0]-irfm['IRFM TEFF'][cl],dr13['FPARAM'][cl,3],zr=[-2.5,0.5],xr=[6000,3000],yr=[-300,300],xt='Teff')
    mean=bindata(dr13['FPARAM'][cl,0][0],(dr13['FPARAM'][cl,0]-irfm['IRFM TEFF'][cl])[0],bins)
    plots.plotp(ax[0,1],bins,mean,marker='o',size=40)
    ax[0,1].text(0.1,0.9,'CL',transform=ax[0,1].transAxes)

    plots.plotc(ax[1,0],dr13['FPARAM'][th,0],dr13['FPARAM'][th,0]-irfm['IRFM TEFF'][th],dr13['FPARAM'][th,3],zr=[-2.5,0.5],xr=[6000,3000],yr=[-300,300],xt='Teff',yt='ASPCAP-photometric Teff')
    mean=bindata(dr13['FPARAM'][th,0][0],(dr13['FPARAM'][th,0]-irfm['IRFM TEFF'][th])[0],bins)
    plots.plotp(ax[1,0],bins,mean,marker='o',size=40)
    ax[1,0].text(0.1,0.9,'TH',transform=ax[1,0].transAxes)

    plots.plotc(ax[1,1],dr13['FPARAM'][sfd,0],dr13['FPARAM'][sfd,0]-irfm['IRFM TEFF'][sfd],dr13['FPARAM'][sfd,3],zr=[-2.5,0.5],xr=[6000,3000],yr=[-300,300],xt='Teff')
    mean=bindata(dr13['FPARAM'][sfd,0][0],(dr13['FPARAM'][sfd,0]-irfm['IRFM TEFF'][sfd])[0],bins)
    plots.plotp(ax[1,1],bins,mean,marker='o',size=40)
    ax[1,1].text(0.1,0.9,'SFD',transform=ax[1,1].transAxes)

    fig.savefig('dteff_teff.jpg')


def dr13dr12() :
    '''
    compare dr13 dr12 Teff
    '''

    apload.dr12()
    dr12=apload.allStar()[1].data
    apload.dr13()
    dr13=apload.allStar()[1].data
    i1,i2 = match.match(dr12['APOGEE_ID'],dr13['APOGEE_ID'])
    dr12=dr12[i1]
    dr13=dr13[i2]

    fig,ax=plots.multi(1,2,hspace=0.001,wspace=0.001)
    plots.plotc(ax[0,0],dr13['M_H'],dr13['TEFF']-dr12['TEFF'],dr13['TEFF'])

    plots.plotc(ax[1,0],dr13['TEFF'],dr13['TEFF']-dr12['TEFF'],dr13['M_H'])
