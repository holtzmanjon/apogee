# routines for comparing gravities with asteroseismic sample

from holtz.apogee import apload
from holtz.apogee import apselect
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

def apokasc(allstar,logg='LOGG_SYD_SCALING',apokasc='APOKASC_cat_v3.6.0.fits',raw=True,cal=True) :
    '''
    asteroseismic log g comparisons for input allStar structure
    '''

    # match ASPCAP with APOKASC, and get RC/RGB stars
    apokasc=fits.open(os.environ['IDLWRAP_DIR']+'/data/'+apokasc)[1].data
    i1,i2=match.match(allstar['APOGEE_ID'],apokasc['2MASS_ID'])
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
    notrc=np.where(apokasc['CONS_EVSTATES'][i2] != 'RC')[0]

    # Do a 2D fit for RGB stars
    ax=plots.ax()
    newfit = fit.fit2d(allstar['FPARAM'][i1[notrc],3],allstar['FPARAM'][i1[notrc],1],
        allstar['FPARAM'][i1[notrc],1]-apokasc[logg][i2[notrc]],zr=[-1,0.5],gdrange=[-2,2],xr=[-3,1],yr=[1,4],degree=2,plot=ax)
    pdb.set_trace()

    # set up plots
    fig,ax=plots.multi(2,3,hspace=0.5,wspace=0.001)
    if not raw or not cal :
        fig,tmpax=plots.multi(1,3,hspace=0.5,wspace=0.001)
        ax[:,0] = tmpax

    # diff color-coded by gravity as f([M/H])
    # diff color-coded by [M/H] as f(log g)
    # RGB and RC as f(log g)
    if raw :
        plots.plotc(ax[0,0],allstar['FPARAM'][i1,3],allstar['FPARAM'][i1,1]-apokasc[logg][i2],
           allstar['FPARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-1,1],xt='[M/H]',yt='ASPCAP-seismic log g')
        plots.plotc(ax[1,0],allstar['FPARAM'][i1,1],allstar['FPARAM'][i1,1]-apokasc[logg][i2],
           allstar['FPARAM'][i1,3],xr=[0,5],zr=[-2.5,0.5],yr=[-1,1],zt='[M/H]',yt='ASPCAP-seismic log g',xt='log g')
        plots.plotp(ax[2,0],allstar['FPARAM'][i1[rgb],1],allstar['FPARAM'][i1[rgb],1]-apokasc[logg][i2[rgb]],
           xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.8,'RGB'],color='r')
        plots.plotp(ax[2,0],allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],1]-apokasc[logg][i2[rc]],
           xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='b')
        #plots.plotc(ax[3,0],allstar['FPARAM'][i1[rgb],1],allstar['PARAM'][i1[rgb],1]-allstar['FPARAM'][i1[rgb],1],
        #   allstar['FPARAM'][i1[rgb],3],xr=[0,5],yr=[-1,1],xt='seismic log g',yt='corrected-raw log g',label=[0.1,0.9,'allstar (Kurucz)'],zr=[-2,0.5])

    if cal :
        param=newfit(allstar['FPARAM'][i1[notrc],3],allstar['FPARAM'][i1[notrc],1])
        plots.plotc(ax[0,1],allstar['FPARAM'][i1[notrc],3],allstar['FPARAM'][i1[notrc],1]-param-apokasc[logg][i2[notrc]],
           allstar['FPARAM'][i1[notrc],1],zr=[0,5],xr=[-2.5,0.5],yr=[-1,1],xt='[M/H]',colorbar=True,zt='log g')

        #plots.plotc(ax[0,1],allstar['PARAM'][i1,3],allstar['PARAM'][i1,1]-apokasc[logg][i2],
        #   allstar['PARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-1,1],xt='[M/H]',colorbar=True,zt='log g')
        plots.plotc(ax[1,1],allstar['PARAM'][i1,1],allstar['PARAM'][i1,1]-apokasc[logg][i2],
           allstar['PARAM'][i1,3],xr=[0,5],zr=[-2.5,0.5],yr=[-1,1],zt='[M/H]',colorbar=True,xt='log g')
        plots.plotp(ax[2,1],allstar['PARAM'][i1[rgb],1],allstar['PARAM'][i1[rgb],1]-apokasc[logg][i2[rgb]],
           xr=[0,5],yr=[-1,1],xt='log g',color='r')
        plots.plotp(ax[2,1],allstar['PARAM'][i1[rc],1],allstar['PARAM'][i1[rc],1]-apokasc[logg][i2[rc]],
           xr=[0,5],yr=[-1,1],xt='log g',color='b')
        #plots.plotc(ax[3,1],allstar['FPARAM'][i1[rc],1],allstar['PARAM'][i1[rc],1]-allstar['FPARAM'][i1[rc],1],
        #    allstar['FPARAM'][i1[rc],3],xr=[0,5],yr=[-1,1],xt='seismic log g',zr=[-2,0.5])

    pdb.set_trace()

def clusters(allstar,xr=[-2.75,0.5],yr=[-1.,1.],zr=[3500,5500],apokasc='APOKASC_cat_v3.6.0.fits',firstgen=False) :
    '''
    Compare ASPCAP gravities in clusters to physical gravities
    '''
    fig,ax=plots.multi(1,2,hspace=0.001)

    # put APOKASC underneath
    apokasc=fits.open(os.environ['IDLWRAP_DIR']+'/data/'+apokasc)[1].data
    i1,i2=match.match(allstar['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0],allstar['FPARAM'][i1,3],allstar['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],allstar['FPARAM'][i1,0],zr=zr)
    plots.plotc(ax[1],allstar['PARAM'][i1,3],allstar['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],allstar['PARAM'][i1,0],zr=zr)

    # physical gravities
    clust=apselect.clustdata()
    itext=0
    out=open('clust.txt','w')
    for cluster in ['M92','M15','M53','M2','M13','M3','M5','N2420','M67','N6819','N6791'] :
        i=np.where(clust.name == cluster)
        dist=clust[i].dist*1000.
        mh=clust[i].mh
        mass=clust[i].giant_mass
        ejk=0.452*clust[i].ebv
        ah=1.55*clust[i].ebv
        age=np.log10(clust[i].age*1.e9)
        name=clust[i].name
        ytext=0.85-itext%3*0.15
        if mass > 0 :
            # get cluster members
            j=np.array(apselect.clustmember(allstar,cluster,raw=True,firstgen=firstgen))

            # calculate physical gravities
            lum=10.**(-0.4*(allstar['H'][j]-ah+isochrones.bc(allstar['FPARAM'][j,0],filt='h',agerange=[age-0.05,age+0.05])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
            logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*allstar['FPARAM'][j,0]**4/lum)
            plots.plotc(ax[0],allstar['FPARAM'][j,3]*0+mh,allstar['FPARAM'][j,1]-logg,allstar['FPARAM'][j,0],xr=xr,yr=yr,zr=zr,yt='ASPCAP-physical log g')
            ax[0].text(0.9,0.1,'raw',transform=ax[0].transAxes,ha='right')

            plots.plotp(ax[0],allstar['FPARAM'][j,3]*0+mh,allstar['FPARAM'][j,1]-logg,color='k')
            plots.plotp(ax[0],mh[0],np.median(allstar['FPARAM'][j,1]-logg),size=40,color='r')
            ax[0].text(mh[0],ytext,name[0],ha='center')

            out.write('{:<20s}{:8.3f}{:8.3f}{:8.3f}\n'.format(clust[i].name[0],clust[i].dist[0],clust[i].ebv[0],mass[0]))

            gd=np.where((allstar['PARAM'][j,3]>-9)&(allstar['PARAM'][j,1]>-9))[0]
            axim=plots.plotc(ax[1],allstar['PARAM'][j[gd],3]*0+mh,allstar['PARAM'][j[gd],1]-logg[gd],allstar['PARAM'][j[gd],0],xr=xr,yr=yr,zr=zr,xt='[M/H]',yt='ASPCAP-physical log g')
            ax[1].text(0.9,0.1,'calibrated',transform=ax[1].transAxes,ha='right')

            plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j[gd],1]-logg[gd]),size=40)
            # apply a temperature correction for the physical gravities
            logg_new=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*(allstar['FPARAM'][j,0]-100.*allstar['FPARAM'][j,3])**4/lum)
            plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j,1]-logg_new),size=40,color='b')
            # use a photometric temperature
            logg_phot=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*stars.ghb(allstar['J'][j]-allstar['K'][j]-ejk,allstar['FPARAM'][j,3])**4/lum)
            plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j,1]-logg_phot),size=40,color='g')
            ax[1].text(mh[0],ytext,name[0],ha='center')
            itext+=1

    # Now adding the colorbar
    cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8]) 
    cb = plt.colorbar(axim, cax = cbaxes)  
    out.close()


def dr13dr12() :
    '''
    ASPCAP compared with physical and asteroseismic log g, DR13/DR12/l30i
    '''
    apokasc=fits.open(os.environ['IDLWRAP_DIR']+'/data/'+apokasc)[1].data
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
        plots.plotc(ax[0,0],dr12['FPARAM'][j,3]*0+mh,dr12['FPARAM'][j,1]-logg,dr12['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],yt='ASPCAP-physical log g',label=[0.1,0.9,'DR2 raw'])
        plots.plotc(ax[1,0],dr12['PARAM'][j,3]*0+mh,dr12['PARAM'][j,1]-logg,dr12['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-physical log g',label=[0.1,0.9,'DR12 cal'])

        #DR13
        j=apselect.clustmember(dr13,cluster,raw=True)
        lum=10.**(-0.4*(dr13['H'][j]+isochrones.bc(dr13['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*dr13['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,1],dr13['FPARAM'][j,3]*0+mh,dr13['FPARAM'][j,1]-logg,dr13['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'DR13 raw'])
        plots.plotc(ax[1,1],dr13['PARAM'][j,3]*0+mh,dr13['PARAM'][j,1]-logg,dr13['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'DR13 cal'],xt='[M/H]')

        #l30i
        j=apselect.clustmember(l30i,cluster,raw=True)
        lum=10.**(-0.4*(l30i['H'][j]+isochrones.bc(l30i['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*l30i['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,2],l30i['FPARAM'][j,3]*0+mh,l30i['FPARAM'][j,1]-logg,l30i['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'l30i raw'])
        plots.plotc(ax[1,2],l30i['PARAM'][j,3]*0+mh,l30i['PARAM'][j,1]-logg,l30i['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'l30i cal'],xt='[M/H]')
    
    plt.show()
    pdb.set_trace()

    # plots vs asterseismic
    fig,ax =plots.multi(3,2,wspace=0.001,hspace=0.001)
    
    i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,0],dr13['FPARAM'][i1,3],dr13['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr13['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[0.1,0.9,'DR13 raw'])
    plots.plotc(ax[1,0],dr13['PARAM'][i1,3],dr13['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr13['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-seismic log g',label=[0.1,0.9,'DR13 cal'])

    i1,i2=match.match(dr12['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,1],dr12['FPARAM'][i1,3],dr12['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr12['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'DR12 raw'])
    plots.plotc(ax[1,1],dr12['PARAM'][i1,3],dr12['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr12['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[0.1,0.9,'DR12 cal'])

    i1,i2=match.match(l30i['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,2],l30i['FPARAM'][i1,3],l30i['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],l30i['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'l30i raw'])
    plots.plotc(ax[1,2],l30i['PARAM'][i1,3],l30i['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],l30i['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[0.1,0.9,'l30i cal'])

    plt.show()

def kurucz_marcs(logg='LOGG_SYD_SCALING',apokasc='APOKASC_cat_v3.6.0.fits') :
    '''
    asteroseismic log g comparisons for Kurucz and MARCS results
    '''
    # read APOKASC
    apokasc=fits.open('APOKASC_cat_v3.6.0.fits')[1].data
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
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
    #plt.plot(apokasc[logg][i2],l30i['FPARAM'][i1,1],'ro')
    # color code by [M/H]
    #plots.plotc(ax1,l30i['FPARAM'][i1,0],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,3],zr=[-1,0.5],xr=[3500,5000],yr=[-1,1],xt='Teff',yt='ASPCAP-seismic log g',label=[0.1,0.9,'l30i (MARCS)'],colorbar=True,zt='[M/H]')
    # color code by [alpha/M]
    plots.plotc(ax1,l30i['FPARAM'][i1,0],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,6],zr=[-0.1,0.4],xr=[3500,5000],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.1,0.9,'l30i (MARCS)'],colorbar=True,zt='[alpha/M]')
    # plot ATMOS_HOLE_WARN and BAD points 
    plots.plotp(ax1,l30i['FPARAM'][i1[warn],0],l30i['FPARAM'][i1[warn],1]-apokasc[logg][i2[warn]],color='y')
    plots.plotp(ax1,l30i['FPARAM'][i1[bad],0],l30i['FPARAM'][i1[bad],1]-apokasc[logg][i2[bad]],color='r')
    # colod code by Teff
    #plots.plotc(ax2,l30i['FPARAM'][i1,1],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,0],zr=[3500,5000],xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.1,0.9,'l30i'],colorbar=True,zt='Teff')

    # match dr13 with APOKASC
    i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
    #plt.plot(apokasc[logg][i2],dr13['FPARAM'][i1,1],'bo')
    #plt.xlim(1.60,1.15)
    #plt.ylim(1.15,1.85)
    #plots.plotc(ax2,dr13['FPARAM'][i1,0],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,3],zr=[-1,0.5],xr=[3500,5000],yr=[-1,1],xt='Teff',yt='ASPCAP-seismic log g',label=[0.1,0.9,'dr13 (Kurucz)'],colorbar=True,zt='[M/H]')
    plots.plotc(ax2,dr13['FPARAM'][i1,0],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,6],zr=[-0.1,0.4],xr=[3500,5000],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.1,0.9,'dr13 (Kurucz)'],colorbar=True,zt='[alpha/M]')
    ##plots.plotc(ax4,dr13['FPARAM'][i1,1],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,0],zr=[3500,5000],xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.1,0.9,'dr13'],colorbar=True,zt='Teff')

    plt.tight_layout()
    plt.show()

