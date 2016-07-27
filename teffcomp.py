# routines for comparing gravities with photometric` sample

from sdss.apogee import apload
from sdss.apogee import apselect
from astropy.io import fits
from holtz.gal import isochrones
from holtz.gal import stars
from holtz.tools import match
from holtz.tools import plots
from holtz.tools import fit
from holtz.apogee import flag
import pdb
import matplotlib.pyplot as plt
import numpy as np
import os
import astropy

def bindata(xdata,ydata,bins) :
    """
    Given xdata, ydata, and bins in x, returns mean of ydata in each of the bins
    """
    mean=bins*0.
    for i in range(len(bins)-1) :
      j=np.where((xdata>bins[i]) & (xdata<bins[i+1]))[0]
      mean[i]=ydata[j].mean() 
    return mean

def ghb(allstar,glatmin=30.,ebvmax=0.03,dwarf=False,trange=[4000,5000],mhrange=[-2.5,0.75]) :
    """
    Compares allstar ASPCPAP Teff with photometric Teff from GHB for sample of stars with GLAT>glatmin and SFD_EBV<ebvmax,
    does fits

    Args:
        allstar   : allStar structure

    Keyword args:
        glatmin (float) : minimum GLAT for sample (default=30)
        ebvmax (float)  : maximum SFD_EBV for sample (default=0.03)
        dwarf (bool)    : use dwarfs and dwarf GHB  (default = False)
    """

    # select data to use
    if dwarf :
        gd=apselect.select(allstar,badval=['STAR_BAD'],teff=trange,mh=mhrange,logg=[3.8,5.0],raw=True)
    else :
        gd=apselect.select(allstar,badval=['STAR_BAD'],teff=trange,mh=mhrange,logg=[0,3.8],raw=True)
    allstar=allstar[gd]
    j=np.where((allstar['GLAT']>glatmin)&(allstar['SFD_EBV']<ebvmax))[0]
    allstar=allstar[j]

    # plot Teff difference against metallicity
    fig,ax=plots.multi(1,1,hspace=0.001,wspace=0.001)
    xr=[-3.0,1.0]
    yr=[-400,300]
    zr=[3500,6000]
    bins=np.arange(-2.5,0.75,0.25)
    # diff color-coded by gravity as f([M/H])
    ghb=stars.ghb(allstar['J']-allstar['K'],allstar['FPARAM'][:,3],dwarf=dwarf)
    plots.plotc(ax,allstar['FPARAM'][:,3],allstar['FPARAM'][:,0]-ghb,allstar['FPARAM'][:,0],zr=zr,xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][:,3],allstar['FPARAM'][:,0]-ghb,bins)
    plots.plotp(ax,bins,mean,marker='o',size=40)
    ax.text(0.1,0.9,'EBV<0.02',transform=ax.transAxes)
    # 1D quadratic fit as a function of metallicity
    x=np.linspace(-3,1,200)
    pfit = fit.fit1d(allstar['FPARAM'][:,3],allstar['FPARAM'][:,0]-ghb,ydata=allstar['FPARAM'][:,0],degree=2)
    plots.plotl(ax,x,pfit(x))
    print pfit
   
    # do some test 2D and 1D fits and plots 
    fig,ax=plots.multi(2,2,hspace=0.5,wspace=0.001)
    ax[0,1].xaxis.set_visible(False)
    ax[0,1].yaxis.set_visible(False)
    pfit = fit.fit2d(allstar['FPARAM'][:,3],allstar['FPARAM'][:,0],allstar['FPARAM'][:,0]-ghb,plot=ax[0,0],zr=[-500,200],xt='[M/H]',yt=['Teff'],zt='$\Delta Teff$')
    print pfit
    pfit = fit.fit1d(allstar['FPARAM'][:,3],allstar['FPARAM'][:,0]-ghb,ydata=allstar['FPARAM'][:,0],plot=ax[1,0],zr=[-500,200],xt='[M/H]',yt='$\Delta Teff$',xr=[-2.7,0.9],yr=[3500,5000])
    print pfit
    pfit = fit.fit1d(allstar['FPARAM'][:,0],allstar['FPARAM'][:,0]-ghb,ydata=allstar['FPARAM'][:,3],plot=ax[1,1],zr=[-500,200],xt='Teff',xr=[3900,5100],yr=[-2.5,0.5])
    print pfit


def irfm(allstar,trange=[4000,5000],mhrange=[-2.5,0.75]) :
    '''
    Compares allstar ASPCPAP Teff with various photometric Teff from JAJ compilation (SAGA, CL, TH, SFD)
    Does fits 

    Args:
        allstar   : allStar structure

    '''

    # select stars
    gd=apselect.select(allstar,badval=['STAR_BAD'],teff=trange,mh=mhrange,raw=True)
    allstar=allstar[gd]

    # get IRFM data
    irfm=fits.open(os.environ['IDLWRAP_DIR']+'/data/irfm_temp.fits')[1].data

    # get the subsamples and match. Note that we have to do this separately for each subsample because some
    #   stars appear in more than one subsample
    saga=np.where(irfm['SOURCE'] == 'SAGA')[0]
    saga1,saga2=match.match(np.chararray.strip(allstar['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID'][saga]))
    cl=np.where(irfm['SOURCE'] == 'CL')[0]
    cl1,cl2=match.match(np.chararray.strip(allstar['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID'][cl]))
    th=np.where(irfm['SOURCE'] == 'TH')[0]
    th1,th2=match.match(np.chararray.strip(allstar['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID'][th]))
    sfd=np.where(irfm['SOURCE'] == 'SFD')[0]
    sfd1,sfd2=match.match(np.chararray.strip(allstar['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID'][sfd]))

    # plot diff color-coded by gravity as f([M/H])
    fig,ax=plots.multi(2,2,hspace=0.001,wspace=0.001)
    xr=[-3.0,1.0]
    yr=[-400,300]
    zr=[3500,6000]
    bins=np.arange(-2.5,0.75,0.25)

    # SAGA
    plots.plotc(ax[0,0],allstar['FPARAM'][saga1,3],allstar['FPARAM'][saga1,0]-irfm['IRFM TEFF'][saga[saga2]],allstar['FPARAM'][saga1,0],zr=zr,xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][saga1,3],allstar['FPARAM'][saga1,0]-irfm['IRFM TEFF'][saga[saga2]],bins)
    plots.plotp(ax[0,0],bins,mean,marker='o',size=40)
    ax[0,0].text(0.1,0.9,'SAGA',transform=ax[0,0].transAxes)

    # CL
    plots.plotc(ax[0,1],allstar['FPARAM'][cl1,3],allstar['FPARAM'][cl1,0]-irfm['IRFM TEFF'][cl[cl2]],allstar['FPARAM'][cl1,0],zr=zr,xr=xr,yr=yr,xt='[M/H]')
    mean=bindata(allstar['FPARAM'][cl1,3],(allstar['FPARAM'][cl1,0]-irfm['IRFM TEFF'][cl[cl2]]),bins)
    plots.plotp(ax[0,1],bins,mean,marker='o',size=40)
    ax[0,1].text(0.1,0.9,'CL',transform=ax[0,1].transAxes)

    # TH
    plots.plotc(ax[1,0],allstar['FPARAM'][th1,3],allstar['FPARAM'][th1,0]-irfm['IRFM TEFF'][th[th2]],allstar['FPARAM'][th1,0],zr=zr,xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][th1,3],(allstar['FPARAM'][th1,0]-irfm['IRFM TEFF'][th[th2]]),bins)
    plots.plotp(ax[1,0],bins,mean,marker='o',size=40)
    ax[1,0].text(0.1,0.9,'TH',transform=ax[1,0].transAxes)

    # SFD
    plots.plotc(ax[1,1],allstar['FPARAM'][sfd1,3],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],allstar['FPARAM'][sfd1,0],zr=zr,xr=xr,yr=yr,xt='[M/H]')
    mean=bindata(allstar['FPARAM'][sfd1,3],(allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]]),bins)
    plots.plotp(ax[1,1],bins,mean,marker='o',size=40)
    ax[1,1].text(0.1,0.9,'SFD',transform=ax[1,1].transAxes)

    fig.savefig('dteff_mh.jpg')

    # plot diff color-coded by gravity as f([M/H])
    fig,ax=plots.multi(2,2,hspace=0.001,wspace=0.001)
    zr=[-2.0,0.5]
    yr=[-400,300]
    xr=[6000,3500]
    bins=np.arange(3500,5500,250)

    # SAGA
    plots.plotc(ax[0,0],allstar['FPARAM'][saga1,0],allstar['FPARAM'][saga1,0]-irfm['IRFM TEFF'][saga[saga2]],allstar['FPARAM'][saga1,3],zr=zr,xr=xr,yr=yr,xt='Teff',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][saga1,0],(allstar['FPARAM'][saga1,0]-irfm['IRFM TEFF'][saga[saga2]]),bins)
    plots.plotp(ax[0,0],bins,mean,marker='o',size=40)
    ax[0,0].text(0.1,0.9,'SAGA',transform=ax[0,0].transAxes)

    # CL
    plots.plotc(ax[0,1],allstar['FPARAM'][cl1,0],allstar['FPARAM'][cl1,0]-irfm['IRFM TEFF'][cl[cl2]],allstar['FPARAM'][cl1,3],zr=zr,xr=xr,yr=yr,xt='Teff')
    mean=bindata(allstar['FPARAM'][cl1,0],(allstar['FPARAM'][cl1,0]-irfm['IRFM TEFF'][cl[cl2]]),bins)
    plots.plotp(ax[0,1],bins,mean,marker='o',size=40)
    ax[0,1].text(0.1,0.9,'CL',transform=ax[0,1].transAxes)

    # TH
    plots.plotc(ax[1,0],allstar['FPARAM'][th1,0],allstar['FPARAM'][th1,0]-irfm['IRFM TEFF'][th[th2]],allstar['FPARAM'][th1,3],zr=zr,xr=xr,yr=yr,xt='Teff',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][th1,0],(allstar['FPARAM'][th1,0]-irfm['IRFM TEFF'][th[th2]]),bins)
    plots.plotp(ax[1,0],bins,mean,marker='o',size=40)
    ax[1,0].text(0.1,0.9,'TH',transform=ax[1,0].transAxes)

    # SFD
    plots.plotc(ax[1,1],allstar['FPARAM'][sfd1,0],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],allstar['FPARAM'][sfd1,3],zr=zr,xr=xr,yr=yr,xt='Teff')
    mean=bindata(allstar['FPARAM'][sfd1,0],(allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]]),bins)
    plots.plotp(ax[1,1],bins,mean,marker='o',size=40)
    ax[1,1].text(0.1,0.9,'SFD',transform=ax[1,1].transAxes)

    fig.savefig('dteff_teff.jpg')

    # do 2D fits with Teff and [M/H], and 1D fits with each

    fig,ax=plots.multi(2,2,hspace=0.5,wspace=0.001)
    ax[0,1].xaxis.set_visible(False)
    ax[0,1].yaxis.set_visible(False)
    pfit = fit.fit2d(ax[0,0],allstar['FPARAM'][sfd1,3],allstar['FPARAM'][sfd1,0],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],plot=True,zr=[-500,200],xt='[M/H]',yt=['Teff'],zt='$\Delta Teff$')
    pfit = fit.fit1d(ax[1,0],allstar['FPARAM'][sfd1,3],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],ydata=allstar['FPARAM'][sfd1,0],plot=True,zr=[-500,200],xt='[M/H]',yt='$\Delta Teff$',xr=[-2.7,0.9],yr=[3500,5000])
    pfit = fit.fit1d(ax[1,1],allstar['FPARAM'][sfd1,0],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],ydata=allstar['FPARAM'][sfd1,3],plot=True,zr=[-500,200],xt='Teff',xr=[3900,5100],yr=[-2.5,0.5])

    pdb.set_trace()

    return pfit


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
    plots.plotc(ax[0],dr13['M_H'],dr13['TEFF']-dr12['TEFF'],dr13['TEFF'],xr=[-2.5,0.75],yr=[-300,300],zr=[3500,5000])

    plots.plotc(ax[1],dr13['TEFF'],dr13['TEFF']-dr12['TEFF'],dr13['M_H'],xr=[6500,3000],yr=[-300,300],zr=[-2,0.5])

