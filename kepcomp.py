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

def kurucz_marcs(logg='LOGG_SYD_SCALING',apokasc='APOKASC_cat_v3.6.0.fits') :
    '''
    asteroseismic log g comparisons for Kurucz and MARCS results
    '''
    # read APOKASC
    apokasc=fits.open(apokasc)[1].data
    #j=np.where((apokasc['TEFF_FIT'] < 4000) & (apokasc[logg] > -500))[0]
    #j=np.where((apokasc['CONS_EVSTATES'] == 'RGB'))[0]
    #apokasc=apokasc[j]

    # read DR13 and l30i
    apload.dr13()
    dr13=apload.allStar()[1].data
    apload.aspcap='l30i'
    apload.results='l30i'
    l30i=apload.allStar()[1].data

    # define axes
    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)
    #ax3=fig.add_subplot(223)
    #ax4=fig.add_subplot(224)

    # match l30i with APOKASC
    i1,i2=match.match(l30i['APOGEE_ID'],apokasc['2MASS_ID'])
    warn=np.where(l30i['ASPCAPFLAG'][i1] & flag.aspcapflagval('ATMOS_HOLE_WARN'))[0]
    bad=np.where(l30i['ASPCAPFLAG'][i1] & flag.aspcapflagval('ATMOS_HOLE_BAD'))[0]
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB'))[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC'))[0]
    #plt.plot(apokasc[logg][i2],l30i['FPARAM'][i1,1],'ro')
    # color code by [M/H]
    #plots.plotc(ax1,l30i['FPARAM'][i1,0],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,3],zr=[-1,0.5],xr=[3500,5000],yr=[-1,1],xt='Teff',yt='ASPCAP-Kepler log g',label=[3600,0.7,'l30i (MARCS)'],colorbar=True,zt='[M/H]')
    # color code by [alpha/M]
    plots.plotc(ax1,l30i['FPARAM'][i1,0],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,6],zr=[-0.1,0.4],xr=[3500,5000],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'l30i (MARCS)'],colorbar=True,zt='[alpha/M]')
    # plot ATMOS_HOLE_WARN and BAD points 
    plots.plotp(ax1,l30i['FPARAM'][i1[warn],0],l30i['FPARAM'][i1[warn],1]-apokasc[logg][i2[warn]],color='y')
    plots.plotp(ax1,l30i['FPARAM'][i1[bad],0],l30i['FPARAM'][i1[bad],1]-apokasc[logg][i2[bad]],color='r')
    # colod code by Teff
    #plots.plotc(ax2,l30i['FPARAM'][i1,1],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,0],zr=[3500,5000],xr=[0,5],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[4,0.7,'l30i'],colorbar=True,zt='Teff')

    # match dr13 with APOKASC
    i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
    #plt.plot(apokasc[logg][i2],dr13['FPARAM'][i1,1],'bo')
    #plt.xlim(1.60,1.15)
    #plt.ylim(1.15,1.85)
    #plots.plotc(ax2,dr13['FPARAM'][i1,0],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,3],zr=[-1,0.5],xr=[3500,5000],yr=[-1,1],xt='Teff',yt='ASPCAP-Kepler log g',label=[3600,0.7,'dr13 (Kurucz)'],colorbar=True,zt='[M/H]')
    plots.plotc(ax2,dr13['FPARAM'][i1,0],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,6],zr=[-0.1,0.4],xr=[3500,5000],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'dr13 (Kurucz)'],colorbar=True,zt='[alpha/M]')
    ##plots.plotc(ax4,dr13['FPARAM'][i1,1],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,0],zr=[3500,5000],xr=[0,5],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[4,0.7,'dr13'],colorbar=True,zt='Teff')

    plt.tight_layout()
    plt.show()

def dr13(logg='LOGG_SYD_SCALING',apokasc='APOKASC_cat_v3.6.0.fits') :
    '''
    asteroseismic log g comparisons for DR13
    '''

    apload.dr13()
    dr13=apload.allStar()[1].data
    apokasc=fits.open(apokasc)[1].data
    i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB'))[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC'))[0]

    fig=plt.figure()
    # diff color-coded by gravity as f([M/H])
    ax=fig.add_subplot(321)
    plots.plotc(ax2,dr13['FPARAM'][i1,3],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-1,1],xt='[M/H]',yt='ASPCAP-Kepler log g',colorbar=True,zt='log g')
    ax=fig.add_subplot(322)
    plots.plotc(ax2,dr13['PARAM'][i1,3],dr13['PARAM'][i1,1]-apokasc[logg][i2],dr13['PARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-1,1],xt='[M/H]',yt='ASPCAP-Kepler log g',colorbar=True,zt='log g')

    # diff color-coded by [M/H] as f(log g)
    ax=fig.add_subplot(323)
    plots.plotc(ax2,dr13['FPARAM'][i1,1],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,3],xr=[0,5],zr=[-2.5,0.5],yr=[-1,1],zt='[M/H]',yt='ASPCAP-Kepler log g',colorbar=True,xt='log g')
    ax=fig.add_subplot(324)
    plots.plotc(ax2,dr13['PARAM'][i1,1],dr13['PARAM'][i1,1]-apokasc[logg][i2],dr13['PARAM'][i1,3],xr=[0,5],zr=[-2.5,0.5],yr=[-1,1],zt='[M/H]',yt='ASPCAP-Kepler log g',colorbar=True,xt='log g')

    # RGB and RC a f(log g)
    ax=fig.add_subplot(325)
    plots.plotp(ax2,dr13['FPARAM'][i1[rgb],1],dr13['FPARAM'][i1[rgb],1]-apokasc[logg][i2[rgb]],xr=[0,5,yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'dr13 (Kurucz)'],colorbar=True,color='r')
    plots.plotp(ax2,dr13['FPARAM'][i1[rc],1],dr13['FPARAM'][i1[rc],1]-apokasc[logg][i2[rc]],xr=[0,5],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'dr13 (Kurucz)'],colorbar=True,color='b')
    ax=fig.add_subplot(326)
    plots.plotp(ax2,dr13['PARAM'][i1[rgb],1],dr13['PARAM'][i1[rgb],1]-apokasc[logg][i2[rgb]],xr=[0,5],yr=[-1,1],xt='log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'dr13 (Kurucz)'],colorbar=True,color='r')
    plots.plotp(ax2,dr13['PARAM'][i1[rc],1],dr13['PARAM'][i1[rc],1]-apokasc[logg][i2[rc]],xr=[0,5],yr=[-1,1],xt='log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'dr13 (Kurucz)'],colorbar=True,color='b')
    


def dr13dr12i() :
    '''
    ASPCAP compared with physical and asteroseismic log g, DR13/DR12/l30i
    '''
    apokasc=fits.open('APOKASC_cat_v3.6.0.fits')[1].data
    j=np.where(apokasc['LOGG_SYD_SCALING'] > -1)[0]
    apokasc=apokasc[j]

    apload.dr12()
    dr12=apload.allStar()[1].data
    apload.dr13()
    dr13=apload.allStar()[1].data
    apload.aspcap='l30i'
    apload.results='l30i'
    l30i=apload.allStar()[1].data

    fig,ax =plots.multi(3,2,wspace=0.001,hspace=0.001)

    # physical gravities
    clust=apselect.clustdata()
    for cluster in ['M92','M15','M53','M2','M13','M3','M5'] :
        i=np.where(clust.name == cluster)
        dist=clust[i].dist*1000.
        mh=clust[i].mh
        mass=0.85

        #DR12
        j=apselect.clustmember(dr12,cluster,raw=True)
        lum=10.**(-0.4*(dr12['H'][j]+isochrones.bc(dr12['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*dr12['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,0],dr12['FPARAM'][j,3]*0+mh,dr12['FPARAM'][j,1]-logg,dr12['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],yt='ASPCAP-physical log g',label=[-2.5,1.5,'DR2 raw'])
        plots.plotc(ax[1,0],dr12['PARAM'][j,3]*0+mh,dr12['PARAM'][j,1]-logg,dr12['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-physical log g',label=[-2.5,1.5,'DR12 cal'])

        #DR13
        j=apselect.clustmember(dr13,cluster,raw=True)
        lum=10.**(-0.4*(dr13['H'][j]+isochrones.bc(dr13['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*dr13['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,1[,dr13['FPARAM'][j,3]*0+mh,dr13['FPARAM'][j,1]-logg,dr13['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[-2.5,1.5,'DR13 raw'])
        plots.plotc(ax[1,1],dr13['PARAM'][j,3]*0+mh,dr13['PARAM'][j,1]-logg,dr13['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[-2.5,1.5,'DR13 cal'],xt='[M/H]')

        #l30i
        j=apselect.clustmember(l30i,cluster,raw=True)
        lum=10.**(-0.4*(l30i['H'][j]+isochrones.bc(l30i['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*l30i['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,2],l30i['FPARAM'][j,3]*0+mh,l30i['FPARAM'][j,1]-logg,l30i['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[-2.5,1.5,'l30i raw'])
        plots.plotc(ax[1,2],l30i['PARAM'][j,3]*0+mh,l30i['PARAM'][j,1]-logg,l30i['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[-2.5,1.5,'l30i cal'],xt='[M/H]')
    
    plt.show()
    pdb.set_trace()

    # plots vs asterseismic
    fig,ax =plots.multi(3,2,wspace=0.001,hspace=0.001)
    
    i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,0],dr13['FPARAM'][i1,3],dr13['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr13['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[-2.,1.5,'DR13 raw'])
    plots.plotc(ax[1,0],dr13['PARAM'][i1,3],dr13['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr13['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-Kepler log g',label=[-2.,1.5,'DR13 cal'])

    i1,i2=match.match(dr12['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,1],dr12['FPARAM'][i1,3],dr12['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr12['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],label=[-2,1.5,'DR12 raw'])
    plots.plotc(ax[1,1],dr12['PARAM'][i1,3],dr12['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr12['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[-2,1.5,'DR12 cal'])

    i1,i2=match.match(l30i['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,2],l30i['FPARAM'][i1,3],l30i['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],l30i['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],label=[-2,1.5,'l30i raw'])
    plots.plotc(ax[1,2],l30i['PARAM'][i1,3],l30i['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],l30i['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[-2,1.5,'l30i cal'])

    plt.show()
