# routines related to APOGEE/ASPCAP calibration

import astropy.io.fits as fits
import os
import numpy as np
import pdb
import matplotlib.pyplot as plt
import plots
import html
import flag
import starflag
import glob
import sdss.apogee.apload as apload
import sdss.apogee.applot as applot

def kepmatch(cal) :
    """
    Given calibration structure, match by ID to Kepler catalog
    Args :
           cal : calibration structure (aspcapField)
    Returns:
           calindex :  list of indices in cal that have Kepler matches
           kep :       Kepler structure
           kepindex :  list of indices in kep that have APOGEE matches
 
    """
    #kep=fits.open(os.getenv('IDLWRAP_DIR')+'/data/APOKASC_Catalog.v7.3.fits')[1].data
    kep=fits.open(os.getenv('IDLWRAP_DIR')+'/data/APOKASC_cat_v3.0.4.fits')[1].data

    from esutil import htm
    h=htm.HTM()
    maxrad=1./3600.
    m1,m2,rad=h.match(cal['RA'],cal['DEC'],kep['RA'],kep['DEC'],maxrad)
    return m1, kep, m2

    #calindex=[]
    #kepindex=[]
    #istar=0
    #for star in kep['2MASS_ID'] :
    #    tmp=star[1:]
    #    j=np.where(cal['APOGEE_ID'].find(tmp) >= 0)[0]
    #    if len(j) > 1 :
    #        print('ERROR, more than one match!')
    #    elif len(j) > 0 :
    #        kepindex.append(istar)
    #        calindex.append(j)
    #    istar+=1
    #return calindex, kep, kepindex



def plot_logg(cal,kep,fit1d=None,fit2d=None, rcfit2d=None) :
    """
    Make lots of log g plots
    """
    logg = cal.FPARAM[:,1]
    teff = cal.FPARAM[:,0]
    mh = cal.FPARAM[:,3]
    #meanfib = cal.MEANFIB
    #diff = logg - kep.KASC_RG_LOGG_SCALE_1
    diff = logg - kep['LOGGRG']
    gd = np.where(kep['LOGGRG'] > 0.1)
    gd=gd[0]

    rgb = np.where(kep['STELLO_EVSTATES'].find('RGB') >= 0)
    rc = np.where(kep['STELLO_EVSTATES'].find('CLUMP') >= 0)
    crgb=np.where(dt(teff[gd],logg[gd],mh[gd]) < 0)
    crgb=gd[crgb[0]]
    crc=np.where(dt(teff[gd],logg[gd],mh[gd]) > 100)
    crc=gd[crc[0]]

    # Figure setup 
    size=(6,4)
    figs=[]
    ytitle=[]

    fig = plt.figure(figsize=size)
    name='plots/logg.jpg'
    ytitle.append('color by evolutionary state')
    ax = fig.add_subplot(1,1,1)
    plots.plotc(ax,logg,diff,'k',xr=[1,3.5],yr=[-0.4,0.5],marker='o',xt='log g',yt='$\Delta$ logg')
    #rgb = np.where(kep['SEISMO EVOL'].find('RGB') >= 0)
    #rc = np.where(kep['SEISMO EVOL'].find('CLUMP') >= 0)
    plots.plotc(ax,logg[crgb],diff[crgb],'m',marker='o')
    plots.plotc(ax,logg[crc],diff[crc],'c',marker='o')
    plots.plotc(ax,logg[rgb],diff[rgb],'r',marker='o')
    plots.plotc(ax,logg[rc],diff[rc],'b',marker='o')
    plt.savefig(name)
    figs.append(['../'+name])
    plt.close()

    fig = plt.figure(figsize=size)
    name='plots/logg_rgb.jpg'
    ytitle.append('color by [M/H]')
    ax = fig.add_subplot(1,1,1)
    plots.plotc(ax,logg,diff,'k',xr=[1,3.5],yr=[-0.4,0.5],marker='o',xt='log g',yt='$\Delta$ logg')
    plots.plotc(ax,logg[rgb],diff[rgb],mh[rgb],zr=[-2.5,0.5],marker='o',colorbar=True,zt='[M/H]')
    plots.plotc(ax,logg[crgb],diff[crgb],mh[crgb],marker='o')
    plt.savefig(name)
    figs.append(['../'+name])
    plt.close()

    fig = plt.figure(figsize=size)
    name='plots/mh.jpg'
    ytitle.append('color by evolutionary state (rgb and crg only)')
    ax = fig.add_subplot(1,1,1)
    plots.plotc(ax,mh,diff,'k',xr=[-2.5,1],yr=[-0.4,0.5],marker='o',xt='[M/H]',yt='$\Delta$ logg')
    plots.plotc(ax,mh[crgb],diff[crgb],'m',marker='o')
    plots.plotc(ax,mh[crc],diff[crc],'c',marker='o')
    plots.plotc(ax,mh[rgb],diff[rgb],'r',marker='o')
    plots.plotc(ax,mh[rc],diff[rc],'b',marker='o')
    if fit1d is not None :
        x=np.arange(-2.5,0.5,0.01)
        plt.plot(x,fit1d(x))
    plt.savefig(name)
    figs.append(['../'+name])
    plt.close()
 
    fig = plt.figure(figsize=size)
    name='plots/logg_mh.jpg'
    ytitle.append('color by [M/H]')
    ax = fig.add_subplot(1,1,1)
    plots.plotc(ax,logg,diff,mh,xr=[1,3.5],yr=[-0.4,0.5],zr=[-2,0.5],cmap='jet_r',xt='log g',yt='$\Delta$ logg',colorbar=True,zt='[M/H]')
    plt.savefig(name)
    figs.append(['../'+name])
    plt.close()

    fig = plt.figure(figsize=size)
    name='plots/mh_logg.jpg'
    ytitle.append('color by logg')
    ax = fig.add_subplot(1,1,1)
    plots.plotc(ax,mh,diff,logg,xr=[-2.5,1],yr=[-0.4,0.5],zr=[1,3.5],cmap='jet_r',xt='[M/H]',yt='$\Delta$ logg',colorbar=True,zt='log g')
    plt.savefig(name)
    figs.append(['../'+name])
    plt.close()

    #fig = plt.figure(figsize=size)
    #name='plots/logg_fiber.jpg'
    #ytitle.append('color by mean fiber')
    #ax = fig.add_subplot(1,1,1)
    #plots.plotc(ax,logg,diff,meanfib,xr=[-1,5],yr=[-0.1,0.5],zr=[0,300],xt='log g',yt='$\Delta$ logg',colorbar=True,zt='<Fiber>')
    #plt.savefig(name)
    #figs.append(['../'+name])
    #plt.close()


    fig = plt.figure(figsize=size)
    name='plots/logg_hr.jpg'
    ytitle.append('color by Delta log g')
    ax = fig.add_subplot(1,1,1)
    plots.plotc(ax,teff[gd],logg[gd],diff[gd],xr=[5500,2500],yr=[5,-1.],zr=[0.,0.4],xt='Teff',yt='logg',colorbar=True,zt='$\Delta$ log g')
    plt.savefig(name)
    figs.append(['../'+name])
    plt.close()

    fig = plt.figure(figsize=size)
    name='plots/logg_dt.jpg'
    ytitle.append('color by Delta log g')
    ax = fig.add_subplot(1,1,1)
    plots.plotc(ax,teff[crgb],logg[crgb],'m',marker='o',xr=[5500,2500],yr=[5,-1.],xt='Teff',yt='logg')
    plots.plotc(ax,teff[crc],logg[crc],'c',marker='o')
    plots.plotc(ax,teff[rgb],logg[rgb],'r',marker='o')
    plots.plotc(ax,teff[rc],logg[rc],'b',marker='o')
    plt.savefig(name)
    figs.append(['../'+name])
    plt.close()

    #plots with fits
    y, x = np.mgrid[0:5:200j, -2.5:1.0:200j]
    if fit2d is not None :
        fig = plt.figure(figsize=size)
        name='plots/logg_rgbfit2d.jpg'
        ytitle.append('color by Delta log g, 2d fit')
        ax = fig.add_subplot(1,1,1)
        plots.plotc(ax,mh[crgb],logg[crgb],diff[crgb],xr=[-2.5,1.0],yr=[5,0.],zr=[0.,0.4],xt='[M/H]',yt='logg',colorbar=True,zt='$\Delta$ log g',linewidth=0)
        ax.imshow(fit2d(x,y),extent=[-2.5,1.0,0,5],
                aspect='auto',vmin=0,vmax=0.4, origin='lower')
        plt.savefig(name)
        figs.append(['../'+name])
        plt.close()

    if rcfit2d is not None :
        fig = plt.figure(figsize=size)
        name='plots/logg_rcfit2d.jpg'
        ytitle.append('color by Delta log g, 2d fit')
        ax = fig.add_subplot(1,1,1)
        plots.plotc(ax,mh[crc],logg[crc],diff[crc],xr=[-2.5,1.0],yr=[5,0.],zr=[0.,0.4],xt='[M/H]',yt='logg',colorbar=True,zt='$\Delta$ log g',linewidth=0)
        ax.imshow(rcfit2d(x,y),extent=[-2.5,1.0,0,5],
                aspect='auto',vmin=0,vmax=0.4, origin='lower')
        plt.savefig(name)
        figs.append(['../'+name])
        plt.close()
  
    if fit1d is not None :
        fig = plt.figure(figsize=size)
        name='plots/logg_hrfit1d.jpg'
        ytitle.append('color by Delta log g, 1d fit')
        ax = fig.add_subplot(1,1,1)
        plots.plotc(ax,teff[gd],logg[gd],diff[gd],xr=[5500,2500],yr=[5,0.],zr=[0.,0.4],xt='Teff',yt='logg',colorbar=True,zt='$\Delta$ log g',linewidth=0,size=10)
        ax.imshow(fit1d(y),extent=[2500,5500,-1,5],
                aspect='auto',vmin=0,vmax=0.4, origin='lower')
        plt.savefig(name)
        figs.append(['../'+name])
        plt.close()

    # make the HTML page with the plots
    html.htmltab(figs,file='html/logg.html', ytitle=ytitle,
                 header='Comparison of ASPCAP 7D results (including vmacro/rot) with APOKASC_Catalog_v7.3')

from astropy.modeling import models, fitting

def fit1d(x,data,degree=1,reject=0) :
    """ Do a 1D fit """
    fit_p = fitting.LevMarLSQFitter()
    p_init = models.Polynomial1D(degree=degree)
    pfit = fit_p(p_init, x, data)
    if reject > 0 :
        gd=np.where(abs(data-pfit(x)) < reject)[0]
        bd=np.where(abs(data-pfit(x)) >= reject)[0]
        #for i in bd :
        #    print(x[i],10.**data[i])
        print('rejected ',len(x)-len(gd),' of ',len(x),' points')
        pfit = fit_p(p_init, x[gd], data[gd])
    return pfit

def fit2d(x,y,data,degree=1,reject=0) :
    """ Do a 2D fit """
    fit_p = fitting.LevMarLSQFitter()
    p_init = models.Polynomial2D(degree=degree)
    pfit = fit_p(p_init, x, y, data)
    if reject > 0 :
        gd=np.where(abs(data-pfit(x,y)) < reject)[0]
        bd=np.where(abs(data-pfit(x,y)) >= reject)[0]
        print( 'rejected ',len(x)-len(gd),' of ',len(x),' points')
        pfit = fit_p(p_init, x[gd], y[gd], data[gd])
    return pfit

def dt(teff,logg,mh) :
    return teff - (4468 + (logg - 2.5)/0.0018 - 382.5*mh )

def fit_logg(cal,kep) :
    """
    Fits for delta log g
    """
    logg = cal.FPARAM[:,1]
    teff = cal.FPARAM[:,0]
    mh = cal.FPARAM[:,3]
    rgb = np.where(kep['STELLO_EVSTATES'].find('RGB') >= 0)
    #meanfib = cal.MEANFIB
    #diff = logg - kep.KASC_RG_LOGG_SCALE_1
    diff = logg - kep.LOGGRG
    gd = np.where((mh>-1) & (logg < 3.8) & (teff<5500) & (teff>3700) & (abs(diff)<0.5) & (kep.LOGGRG > 0.1))[0]
    fit = np.where((kep[gd].STELLO_EVSTATES.find('RGB') >= 0) | (dt(teff[gd],logg[gd],mh[gd]) < 0) )[0]
    print( kep[gd[fit]].STELLO_EVSTATES)
    gfit1d = fit1d(mh[gd[fit]], diff[gd[fit]])
    gfit2d = fit2d(mh[gd[fit]], logg[gd[fit]], diff[gd[fit]])

    fit = np.where((kep[gd].STELLO_EVSTATES.find('RC') >= 0) | (dt(teff[gd],logg[gd],mh[gd]) > 100) )[0]
    print(kep[gd[fit]].STELLO_EVSTATES)
    grcfit1d = fit1d(mh[gd[fit]], diff[gd[fit]])
    grcfit2d = fit2d(mh[gd[fit]], logg[gd[fit]], diff[gd[fit]])
    return gfit1d, gfit2d, grcfit1d, grcfit2d

def fit_vmicro(file,reject=0.) :
    """ Fit microturbulence relation """
    data = fits.open(file)[1].data
    teff, logg, mh, vmicro, vmacro, meanfib, aspcapflag, chi2 = read(data)
    vrange = [0.45,3.5]
    gd = np.where((mh>-1.0) & (logg < 3.8) & (10.**vmicro < vrange[1]) & ((aspcapflag & flag.badaspcapflag()) == 0))
    # gd = np.where((mh>-1.5) & (10.**vmicro < vrange[1]) & ((aspcapflag & flag.badaspcapflag()) == 0))
    vfit1d = fit1d(logg[gd], vmicro[gd],degree=3,reject=reject)
    vfit2d = fit2d(teff[gd], logg[gd], vmicro[gd],degree=3)
    plot_vmicro(teff, logg, mh, meanfib, vmicro, aspcapflag, vrange, vfit1d, vfit2d, vt='vmicro')
    return vfit1d, vfit2d

def fit_vmacro(file,reject=0) :
    """ Fit macroturbulence relation """
    data = fits.open(file)[1].data
    teff, logg, mh, vmicro, vmacro, meanfib, aspcapflag, chi2 = read(data)
    vrange = [0,15]
    # gd = np.where((mh>-1) & (logg < 3.8) & (10.**vmicro < vrange[1]) & ((aspcapflag & flag.badaspcapflag()) == 0))
    gd = np.where((mh>-1.5) & (logg < 3.8) & (10.**vmacro < vrange[1]) & (chi2 <20) & ((aspcapflag & flag.badaspcapflag()) == 0))
    vfit1d = fit1d(logg[gd], vmacro[gd],degree=3,reject=reject)
    #vfit2d = fit2d(np.log10(teff[gd]), logg[gd], vmacro[gd],reject=reject)
    #plot_vmicro(teff, logg, mh, meanfib, vmacro, aspcapflag, vrange, vfit1d, vfit2d, vt='vmacro',logte=True)
    gd = np.where((mh>-1.5) & (logg < 3.8) & (10.**vmacro < vrange[1]) & (chi2 <20) & ((aspcapflag & flag.badaspcapflag()) == 0))
    vfit2d = fit2d(mh[gd], logg[gd], vmacro[gd],reject=reject)
    plot_vmicro(mh, logg, teff, meanfib, vmacro, aspcapflag, vrange, vfit1d, vfit2d, vt='vmacro',xr=[-2.5,1], xt='[M/H]',zr=[2500,5500])
    return vfit1d, vfit2d

def plot_vmicro(teff, logg, mh, meanfib, v, aspcapflag, vr, fit1d, fit2d,
         xr=[5500,2500], yr=[5,-1], zr=[-2.5,0.5],vt='vmicro', logte=False, xt='Teff') :

    gd = np.where( ((aspcapflag & flag.badaspcapflag()) == 0))
    bd = np.where( ((aspcapflag & flag.badaspcapflag()) > 0))

    fig = plt.figure(figsize=(14,8))
    y, x = np.mgrid[yr[1]:yr[0]:200j, xr[1]:xr[0]:200j]
    if logte : 
        x = np.log10(x) 
    ax=fig.add_subplot(3,2,1)
    plots.plotc(ax,teff[gd],logg[gd],10.**v[gd],xr=xr,yr=yr,zr=vr,
                xt=xt,yt='log g',zt=vt,colorbar=True,size=15)
    plots.plotc(ax,teff[bd],logg[bd],10.**v[bd],xr=xr,yr=yr,zr=vr,
                xt=xt,yt='log g',zt=vt,size=15,marker='x',linewidth=1)
    ax.imshow(10.**fit2d(x,y),extent=[xr[1],xr[0],yr[1],yr[0]],
                aspect='auto',vmin=0,vmax=vr[1], origin='lower')

    ax=fig.add_subplot(3,2,2)
    plots.plotc(ax,teff[gd],logg[gd],10.**v[gd],xr=xr,yr=yr,zr=vr,
                xt=xt,yt='log g',zt=vt,colorbar=True,size=15)
    plots.plotc(ax,teff[bd],logg[bd],10.**v[bd],xr=xr,yr=yr,zr=vr,
                xt=xt,yt='log g',zt=vt,size=15,marker='x',linewidth=1)
    ax.imshow(10.**fit1d(y),extent=[xr[1],xr[0],yr[1],yr[0]],
              aspect='auto',vmin=0,vmax=vr[1], origin='lower')

    ax=fig.add_subplot(3,2,3)
    if logte : 
        temp = np.log10(teff) 
    else :
        temp = teff
    plots.plotc(ax, logg, 10.**fit2d(temp,logg), mh, xr=[-1,5], yr=vr, zr=zr,
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')
    #ax.set_yscale('log')
    #ax.set_ylim(vr[0],vr[1])

    ax=fig.add_subplot(3,2,4)
    plots.plotc(ax, logg, 10.**fit1d(logg), mh, xr=[-1,5], yr=vr, zr=zr,
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')
    #ax.set_yscale('log')
    #ax.set_ylim(vr[0],vr[1])

    ax=fig.add_subplot(3,2,5)
    plots.plotc(ax, logg[gd], 10.**v[gd], mh[gd], xr=[-1,5], yr=vr, zr=zr,
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')
    plots.plotc(ax, logg[bd], 10.**v[bd], mh[bd], xr=[-1,5], yr=vr, zr=zr,
                size=10,xt='log g', yt=vt, zt='[M/H]',marker='x',linewidth=1)
    tmp=np.arange(-1,5,0.1)
    ax.plot(tmp, 10.**fit1d(tmp))
    #ax.set_yscale('log')
    #ax.set_ylim(vr[0],vr[1])

    ax=fig.add_subplot(3,2,6)
    plots.plotc(ax, logg[gd], 10.**v[gd], meanfib[gd], xr=[-1,5], yr=vr, zr=[0,300],
                size=10,colorbar=True,xt='log g', yt=vt, zt='<fiber>')
    plots.plotc(ax, logg[bd], 10.**v[bd], meanfib[bd], xr=[-1,5], yr=vr, zr=[0,300],
                size=10,xt='log g', yt=vt, zt='<fiber>',marker='x',linewidth=1)
    #ax.set_yscale('log')
    #ax.set_ylim(vr[0],vr[1])

    plt.show()


def read(data) :
    vmicro = data['fparam'][:,2]
    vmacro = data['fparam'][:,7]
    teff = data['fparam'][:,0]
    logg = data['fparam'][:,1]
    mh = data['fparam'][:,3]
    flag = data['aspcapflag']
    chi2 = data['param_chi2']
    try :
       meanfib = data['meanfib']
    except :
       meanfib = 0.*vmicro
    return teff, logg, mh, vmicro, vmacro, meanfib, flag, chi2


def main() :
    cal=fits.open('cal.fits')[1].data
    calindex, kep, kepindex = kepmatch(cal)
    gfit1d, gfit2d = fit_logg(cal[calindex],kep[kepindex])
    pdb.set_trace()
    plot_logg(cal[calindex],kep[kepindex],fit1d=gfit1d,fit2d=gfit2d)

def plot_hr(file,dochi2=False) :
    global ax

    data = fits.open(file)[1].data
    teff, logg, mh, vmicro, vmacro, meanfib, aspcapflag, chi2 = read(data)
    gd = np.where( ((aspcapflag & flag.badaspcapflag()) == 0))
    bd = np.where( ((aspcapflag & flag.badaspcapflag()) > 0))
    fig = plt.figure(figsize=(14,8))
    plots.event(fig)
    plots._data=data
    plots._data_x=teff
    plots._data_y=logg
    ax=fig.add_subplot(111)
    if dochi2 :
        plots.plotc(ax,teff[bd],logg[bd],chi2[bd],xr=[6500,2500],yr=[5,-1],zr=[0,50],xt='Teff',yt='log g',zt='chi2',size=15,marker='x',linewidth=1)
        plots.plotc(ax,teff[gd],logg[gd],chi2[gd],xr=[6500,2500],yr=[5,-1],zr=[0,50],xt='Teff',yt='log g',zt='chi2',size=15,colorbar=True)
    else :
        plots.plotc(ax,teff[bd],logg[bd],mh[bd],xr=[6500,2500],yr=[5,-1],zr=[-2,0.5],xt='Teff',yt='log g',zt='[M/H]',size=15,marker='x',linewidth=1)
        plots.plotc(ax,teff[gd],logg[gd],mh[gd],xr=[6500,2500],yr=[5,-1],zr=[-2,0.5],xt='Teff',yt='log g',zt='[M/H]',size=15,colorbar=True)
    return ax, data

def plot_hr_fiber(file) :
    data = fits.open(file)[1].data
    teff, logg, mh, vmicro, vmacro, meanfib, aspcapflag, chi2 = read(data)
    gd = np.where( ((aspcapflag & flag.badaspcapflag()) == 0))
    bd = np.where( ((aspcapflag & flag.badaspcapflag()) > 0))
    fig = plt.figure(figsize=(14,8))
    ax=fig.add_subplot(111)
    plots.plotc(ax,teff[bd],logg[bd],meanfib[bd],xr=[6500,2500],yr=[5,-1],zr=[0,300],xt='Teff',yt='log g',zt='<fiber>',size=15,marker='x',linewidth=1)
    plots.plotc(ax,teff[gd],logg[gd],meanfib[gd],xr=[6500,2500],yr=[5,-1],zr=[0,300],xt='Teff',yt='log g',zt='<fiber>',size=15,colorbar=True)

def plot_hr_grid(file) :
    data = fits.open(file)[1].data
    teff, logg, mh, vmicro, vmacro, meanfib, aspcapflag, chi2 = read(data)
    gd = np.where( ((aspcapflag & flag.badaspcapflag()) == 0))
    bd = np.where( ((aspcapflag & flag.badaspcapflag()) > 0))
    gk = np.where( data['CLASS'].find('GK') >= 0)
    m = np.where( data['CLASS'].find('M') >= 0)
    f = np.where( data['CLASS'].find('F') >= 0)
    print(type(teff))
    grid = teff*0.
    grid[f] = 0
    grid[gk] = 1
    grid[m] = 2
    fig = plt.figure(figsize=(14,8))
    ax=fig.add_subplot(111)
    plots.plotc(ax,teff[bd],logg[bd],grid[bd],xr=[6500,2500],yr=[5,-1],zr=[0,2],xt='Teff',yt='log g',zt='grid',size=15,marker='x',linewidth=1)
    plots.plotc(ax,teff[gd],logg[gd],grid[gd],xr=[6500,2500],yr=[5,-1],zr=[0,2],xt='Teff',yt='log g',zt='grid',size=15,colorbar=True)


def plot_spec(apogee_id) :
    # find the directory assuming it is under the current directory
    dir = glob.glob('*/*'+apogee_id+'*')[0].split('/')[0]
    hdu=apload.apStar('../cal/'+dir,apogee_id)
    print(starflag.starflag(hdu[0].header['STARFLAG']))
    hdu=apload.aspcapStar(dir,apogee_id)
    applot.aspcapStar(hdu)
