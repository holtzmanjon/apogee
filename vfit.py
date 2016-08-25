from astropy.io import fits
from astropy.io import ascii
from astropy.modeling import models, fitting
import astropy.constants as const
import matplotlib.pyplot as plt
import numpy as np
import pdb
import matplotlib as mpl
import sys
import argparse
import os
from holtz.tools import fit
from holtz.tools import plots
from holtz.tools import match

def fit_vmicro(file,teffrange=[3550,5500],mhrange=[-2.5,1],loggrange=[-0.5,4.75],vrange=[0,4],vmrange=[0,6],maxerr=0.1, degree=1,reject=0) :
    """ 
    Fit microturbulence relation  with 1D f(log g) and 2D f(Teff, logg) fits, plots

    Args:
        file : file with calib structure

    Keyword args :
        degree : degree of fit (default=1)
        mhrange : range of [M/H] (default=[-1,1])
        loggrange : range of log g (default=[-1,3.8])
        maxerr : maximum uncertainy in vmicro
        vrange  : scaling range for vmacro plot, vrange[1] sets maximum good vmacro

    Returns:
        fit1d, fit2d : 1D and 2D polynomial fits
    """
    data=fits.open(file)[1].data
    vmicro = data['fparam'][:,2]
    vmacro = data['fparam'][:,7]
    # fix locked vmacro by hand (bad!)
    j=np.where(np.isclose(vmacro,1.))[0]
    vmacro[j] = 0.6
    teff = data['fparam'][:,0]
    logg = data['fparam'][:,1]
    mh = data['fparam'][:,3]
    try :
        meanfib = data['meanfib']
    except :
        meanfib = 0.*vmicro
    try :
        ninst= data['NINST'][:,1]-data['NINST'][:,2]
    except :
        ninst= data['NVISITS']*0

    gd = np.where((mh>mhrange[0]) & (mh<mhrange[1]) & (logg > loggrange[0]) & (logg < loggrange[1]) & 
                  (teff>teffrange[0]) & (teff<teffrange[1]) & (10.**vmacro>vmrange[0]) & (10.**vmacro<vmrange[1]) &
                  (np.sqrt(data['fparam_cov'][:,2,2]) < maxerr) & (10.**vmicro > vrange[0]) & (10.**vmicro < vrange[1]))[0]

    # remove non-1st generation GC stars
    gcstars = ascii.read(os.environ['IDLWRAP_DIR']+'/data/gc_szabolcs.dat')
    bd=np.where(gcstars['pop'] != 1)[0]
    gd = [x for x in gd if data[x]['APOGEE_ID'] not in gcstars['id'][bd]]

    # 1D plots a f(log g)
    fig,ax = plots.multi(2,3,figsize=(12,12))
    fit1d = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[0,0],ydata=mh[gd],log=True,xt='log g',yt='vmicro ([M/H])',yr=[-2.5,0.5],colorbar=True,zt='[M/H]')
    # plot ALL points (even outside of fit range)
    plots.plotc(ax[0,0],logg,10.**vmicro,mh,zr=[-2.5,0.5],xr=[-1,5])

    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[0,1],ydata=teff[gd],log=True,xt='log g',yt='vmicro',yr=[3500,5500],pfit=fit1d,colorbar=True,zt='Teff')
    plots.plotc(ax[0,1],logg,10.**vmicro,teff,zr=[3500,5500],xr=[-1,5])
    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[1,0],ydata=meanfib[gd],log=True,xt='log g',yt='vmicro',yr=[0,300],pfit=fit1d,colorbar=True,zt='mean fiber')
    plots.plotc(ax[1,0],logg,10.**vmicro,meanfib,zr=[0,300],xr=[-1,5])
    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[1,1],ydata=10.**vmacro[gd],log=True,xt='log g',yt='vmicro',yr=[0,15],pfit=fit1d,colorbar=True,zt='vmacro')
    plots.plotc(ax[1,1],logg,10.**vmicro,10.**vmacro,zr=[0,15],xr=[-1,5])

    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[2,0],ydata=ninst[gd],log=True,xt='log g',yt='vmicro',yr=[-1,1],pfit=fit1d,colorbar=True,zt='ninst1-ninst2')
    #plots.plotc(ax[3,1],logg[gd1],10.**vmicro[gd1],mh[gd1],zr=[-2.5,0.5],xr=[-1,5])
    # 2d plot
    #junk = fit.fit1d(logg, vmicro,degree=degree,reject=reject,plot=ax[2,0],plot2d=True,ydata=teff,log=True,yt='Teff',xt='log g',xr=[5,-0.5],yr=[6000,3000],pfit=fit1d,zr=[0,4])


    # plot with DR13 relation
    dr13fit=models.Polynomial1D(degree=3)
    dr13fit.parameters=[0.226,-0.0228,0.0297,-0.013]
    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[2,1],ydata=mh[gd],pfit=dr13fit,log=True,xt='log g',yt='vmicro',colorbar=True,zt='[M/H]')

    fig.tight_layout()

    # 2D plots a f(teff, logg)
    #fit2d = fit.fit2d(teff[gd], logg[gd], vmicro[gd],degree=degree,plot=ax[2,0],yr=[5,0],xr=[6000,3500],xt='TTeff',yt='log g')
    #summary plots
    #plot(teff, logg, mh, meanfib, vmicro, vrange, fit1d, fit2d, vt='vmicro')

    return fit1d

def fit_vmacro(file,mhrange=[-1.,1.], loggrange=[-1,3.8], vrange=[1,15],degree=1,apokasc='APOKASC_cat_v3.6.0.fits') :
    """ 
    Fit macroturbulence relation  with 1D f(log g) and 2D f(Teff, logg) fits, plots

    Args:
        file : file with calib structure

    Keyword args :
        degree : degree of fit (default=1)
        mhrange : range of [M/H] (default=[1.,1])
        loggrange : range of log g (default=[1.,3.8])
        vrange  : scaling range for vmacro plot, vrange[1] sets maximum good vmacro

    Returns:
        fit1d, fit2d : 1D and 2D polynomial fits
    """

    data=fits.open(file)[1].data
    vmicro = data['fparam'][:,2]
    teff = data['fparam'][:,0]
    logg = data['fparam'][:,1]
    mh = data['fparam'][:,3]
    try :
       meanfib = data['meanfib']
    except :
       meanfib = 0.*vmicro
    gd = np.where((mh>mhrange[0]) & (mh<mhrange[1]) & (logg > loggrange[0]) (logg < loggrange[1]) & (10.**vmacro < vrange[1]))[0]
    print('len(gd)', len(gd))

    cal=fits.open(file)[1].data
    apokasc=fits.open(os.environ['IDLWRAP_DIR']+'/data/'+apokasc)[1].data
    i1,i2=match.match(cal['APOGEE_ID'],apokasc['2MASS_ID'])
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
    type=mh*0.
    type[i1[rgb]]=1
    type[i1[rc]]=-1

    fig,ax=plots.multi(1,2)
    fit1d = fit.fit1d(logg[gd], vmacro[gd],degree=degree,plot=ax[0],xt='log g', yt='vmacro (mh)',ydata=mh[gd])
    fit1d = fit.fit1d(logg[gd], vmacro[gd],degree=degree,plot=ax[1],xt='log g', yt='vmacro (RGB/RC)',ydata=type[gd])

    fig,ax=plots.multi(1,2)
    fit2d = fit.fit2d(logg[gd], mh[gd], vmacro[gd],degree=degree,plot=ax[0],xt='log g', yt='[M/H]',zt='log(vmacro)')
    # DR13 fit 
    fit2d.parameters=[0.741,-0.0998,-0.225]
    fit2d = fit.fit2d(logg[gd], mh[gd], vmacro[gd],degree=degree,plot=ax[1],xt='log g', yt='[M/H]',zt='log(vmacro)',pfit=fit2d)
    pdb.set_trace()
    plot(teff, logg, mh, meanfib, vmacro, vrange, fit1d, fit2d, vt='vmacro')
    return fit1d, fit2d

def litplot(mass=1, feh=0) : 
    '''
    Plots literature vmacro relations
    '''
#    parser = argparse.ArgumentParser()
#    parser.add_argument("--mass",default=1.)
#    parser.add_argument("--feh",default=0.)
#    args=parser.parse_args()
#    print('Using mass: ', args.mass)
#    mass=float(args.mass)
#    feh=float(args.feh)
#    print('Using feh: ', feh)
#    mass = 0.9

    # setup plots
    mpl.rcParams['xtick.labelsize'] = 8
    mpl.rcParams['ytick.labelsize'] = 8
    p=plt.figure()

    p,ax=plots.multi(2,2)
    lw=0

    # get isochrone data
    for iso in [0.,-1.,-2.] :
        if iso >=0 : 
            c='p' 
        else: 
            c='m'
        cfeh='z'+c+'{:02d}'.format(int(round(abs(iso)*10)))
        print(cfeh)
        a=ascii.read(os.getenv('ISOCHRONE_DIR')+'/'+cfeh+'.dat')
        age=a['col2']
        j=np.where((age == 9.6) | (age == 9.0) | (age == 10.0))
        logl=a['col5'][j]
        iso_te=10.**a['col6'][j]
        iso_logg=a['col7'][j]
        iso_feh=iso_logg*0.+iso

        vm=vmacro_massarotti(iso_te,iso_logg,iso_feh)
        plots.plotc(ax[0,0],iso_te,iso_logg,vm,xr=[7500,2500],yr=[5,-1.],zr=[0,10],linewidth=lw,size=15)
        vm=vmacro_shetrone(iso_te,iso_logg,iso_feh)
        plots.plotc(ax[1,0],iso_te,iso_logg,vm,xr=[7500,2500],yr=[5,-1.],zr=[0,10],linewidth=lw,size=15)
        vm=vmacro_bergemann(iso_te,iso_logg,iso_feh)
        plots.plotc(ax[0,1],iso_te,iso_logg,vm,xr=[7500,2500],yr=[5,-1.],zr=[0,10],linewidth=lw,size=15)
        vm=vmacro_bergemann_new(iso_te,iso_logg,iso_feh)
        plots.plotc(ax[1,1],iso_te,iso_logg,vm,xr=[7500,2500],yr=[5,-1.],zr=[0,10],linewidth=lw,size=15)


    # setup grids for Teff and logg
    o=np.ones([61,51])
    teff=o*(np.arange(51)*100.)+2500
    logg=np.transpose(np.transpose(o)*(np.arange(61)*.1))-1.
    feh=o*0.+feh
    # L = 4*pi*R^2*sigma*Teff^4
    # R^2 = L / (4 pi sigma Teff^4)
    # log R^2 = log L - log (4*pi*sigma) - 4 logte
    # logg = log G + log M - log L + log(4*pi*sigma) + 4 log te
    sigma=const.sigma_sb.cgs.value
    G=const.G.cgs.value
    lsun=const.L_sun.cgs.value
    msun=const.M_sun.cgs.value
    m=mass*msun
    logl=4*np.log10(teff)-logg-np.log10(lsun)+np.log10(4*np.pi*sigma*G)+np.log10(m)

    a=vmacro_massarotti(teff,logl,feh)
    ax[0,0].set_title('Massarotti vmacro')
    axim=ax[0,0].imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500,2500.,5.,-1.],aspect='auto',origin='upper')
    #p.colorbar(axim)

    a=vmacro_shetrone(teff,logl,feh)
    ax[0,1].set_title('Shetrone vmacro')
    axim=ax[0,1].imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500,2500.,5.,-1.],aspect='auto',origin='upper')
    #p.colorbar(axim)
    
    a=vmacro_bergemann(teff,logg,feh)
    ax[1,0].set_title('Bergemann vmacro')
    axim=ax[1,0].imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500.,2500.,5.,-1.],aspect='auto',origin='upper')
    #p.colorbar(axim)

    a=vmacro_bergemann_new(teff,logg,feh)
    ax[1,1].set_title('Bergemann vmacro new')
    axim=ax[1,1].imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500.,2500.,5.,-1.],aspect='auto',origin='upper')
    #p.colorbar(axim)

    p.savefig('vmacro_'+str(mass)+'.jpg')


def vmacro_massarotti(teff,logl,feh) :
    '''
    Returns Massarotti vmacro as f(teff,logl,feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmacro
    '''
    logte=np.log10(teff)
    vm=10**(3.5*logte+0.25*logl-12.67)
    return vm

def vmacro_shetrone(teff,logl,feh) :
    '''
    Returns Shetrone vmacro as f(teff, logl, feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmacro
    '''
    logte=np.log10(teff)
    return 10**(1.85*logte+0.22*logl-6.49)

def vmacro_bergemann(teff,logg,feh) :
    '''
    Returns Bergemann vmacro as f(teff, logl, feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmacro
    '''
    tref=5250. ; to = 5500. ; go = 4.0
    t=teff ; g=logg

    # set up output, which we will fill in pieces for different regimes
    vt=teff*0.

    # t >= Tref and logg >=3.5
    j= np.where((teff>tref) & (logg >= 3.5))
    a1 =  1.15 ; b1 =  7.e-4; c1 = 1.2e-6
    b2 = -0.13 ; c2 =  0.13
    b3 = -0.37 ; c3 = -0.07
    vt[j] = 3*(a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)

    # t < Tref and logg >=3.5
    j= np.where((teff<tref) & (logg >= 3.5))
    a1 =  1.15 ; b1 =  2.e-4; c1 = 3.95e-7
    b2 = -0.13 ; c2 =  0.13
    b3 =  0.0  ; c3 =  0.0
    vt[j] = 3*(a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)

    # logg < 3.5
    j= np.where(logg < 3.5)
    a1 =  1.15 ; b1 =  2.2e-5 ; c1 = -0.5e-7
    b2 = -0.1 ; c2 =  0.04 
    b3 = -0.37  ; c3 = -0.07
    vt[j] = 3*(a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)
    return vt

def vmacro_bergemann_new(teff,logg,feh) :
    '''
    Returns newer Bergemann vmacro as f(teff, logl, feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmacro
    '''

    t0  = 5500
    g0  = 4.0
    vt=teff*0.

    # MS and RGB, Teff >=5250
    j= np.where((teff>5250.) & (logg >= 3.5))
    Rxa = ( 1.15 + 7e-4*(teff[j]-t0) + 1.2e-6*(teff[j]-t0)**2 
            - 0.13*(logg[j]-g0) +     0.13*(logg[j]-g0)**2 
            -    0.37*feh[j]     -    0.07*feh[j]**2 )   
    vt[j] = 3*Rxa

    # MS, Teff <=5250
    j= np.where((teff<5250.) & (logg >= 3.5))
    Rxa = (1.15 + 2e-4*(teff[j]-t0) + 3.95e-7*(teff[j]-t0)**2 
            - 0.13*(logg[j]-g0) +     0.13*(logg[j]-g0)**2 )
    vt[j] = 3*Rxa

    # RGB/AGB
    j= np.where(logg < 3.5)
    Rxa = ( 1.15 + 2.2e-5*(teff[j]-t0) - 0.5e-7*(teff[j]-t0)**2 
            -    0.1*(logg[j]-g0) +   0.04*(logg[j]-g0)**2  
            -    0.37*feh[j]     -    0.07*feh[j]**2 )
    vt[j] = 3*Rxa
    return vt

def vmicro_bergemann(teff,logg,feh) :
    '''
    Returns Bergemann vmicro as f(teff, logl, feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmicro
    '''

    # reference variables
    tref=5500. ; to = 5500. ; go = 4.0
    t=teff ; g=logg

    # set up output variable, which we will fill in pieces
    vt=teff*0.

    # t >= Tref and logg >=3.5
    j= np.where((teff>=tref) & (logg >= 3.5))
    a1 =  1.1 ; b1 =  2.e-4; c1 = 7.95e-7
    b2 = -0.13 ; c2 =  0.13
    b3 = -0.27 ; c3 = -0.1
    vt[j] = (a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)

    # t < Tref and logg >=3.5
    j= np.where((teff<tref) & (logg >= 3.5))
    a1 =  1.1 ; b1 =  2.e-4; c1 = 3.95e-7
    b2 = -0.13 ; c2 =  0.13
    b3 =  -0.6  ; c3 =  -0.2
    vt[j] = (a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)

    # logg < 3.5
    j= np.where(logg < 3.5)
    a1 =  1.3 ; b1 =  2.2e-5 ; c1 = -0.5e-7
    b2 = -0.1 ; c2 =  0.04 
    b3 =  0.0  ; c3 = 0.0
    vt[j] = (a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)
    return vt

def plot(teff, logg, mh, meanfib, v, vr, fit1d, fit2d, 
         xr=[5500,3500], yr=[5,0], vt='vmicro') :
    '''
    auxiliary plotting routine for vmicro, vmacro plots
    '''
    fig = plt.figure(figsize=(14,8))
    y, x = np.mgrid[yr[1]:yr[0]:200j, xr[1]:xr[0]:200j]
    ax=fig.add_subplot(3,2,1)
    plots.plotc(ax,teff,logg,10.**v,xr=xr,yr=yr,zr=vr,
                xt='Teff',yt='log g',zt=vt,colorbar=True,size=15)
    ax.imshow(10.**fit2d(x,y),extent=[xr[1],xr[0],yr[1],yr[0]],
                aspect='auto',vmin=0,vmax=vr[1], origin='lower')

    ax=fig.add_subplot(3,2,2)
    plots.plotc(ax,teff,logg,10.**v,xr=xr,yr=yr,zr=vr,
                xt='Teff',yt='log g',zt=vt,colorbar=True,size=15)
    ax.imshow(10.**fit1d(y),extent=[xr[1],xr[0],yr[1],yr[0]],
              aspect='auto',vmin=0,vmax=vr[1], origin='lower')

    ax=fig.add_subplot(3,2,3)
    plots.plotc(ax, logg, 10.**fit2d(teff,logg), mh, xr=[0,5], yr=vr, zr=[-2.5,0.5],
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')

    ax=fig.add_subplot(3,2,4)
    plots.plotc(ax, logg, 10.**fit1d(logg), mh, xr=[0,5], yr=vr, zr=[-2.5,0.5],
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')
    x=np.arange(0,5,0.01)
    plots.plotl(ax,x, 10.**(0.226 - 0.0228 *x + 0.0297 *x**2 - 0.0113 *x**3))

    ax=fig.add_subplot(3,2,5)
    plots.plotc(ax, logg, 10.**v, mh, xr=[0,5], yr=vr, zr=[-2.5,0.5],
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')

    ax=fig.add_subplot(3,2,6)
    plots.plotc(ax, logg, 10.**v, meanfib, xr=[0,5], yr=vr, zr=[0,300],
                size=10,colorbar=True,xt='log g', yt=vt, zt='<fiber>')

    plt.show()

