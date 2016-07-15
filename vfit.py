from astropy.io import fits
from astropy.io import ascii
import astropy.constants as const
import matplotlib.pyplot as plt
import plots
import numpy as np
import pdb
import matplotlib as mpl
import sys
import argparse
import os

def vmacro_massarotti(teff,logl,feh) :
    logte=np.log10(teff)
    vm=10**(3.5*logte+0.25*logl-12.67)
    return vm

def vmacro_shetrone(teff,logl,feh) :
    logte=np.log10(teff)
    return 10**(1.85*logte+0.22*logl-6.49)

def vmacro_bergemann(teff,logg,feh) :
    tref=5250.
    to = 5500.
    go = 4.0
    vt=teff*0.
    t=teff
    g=logg

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
    tref=5500.
    to = 5500.
    go = 4.0
    vt=teff*0.
    t=teff
    g=logg

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


# fits
from astropy.modeling import models, fitting

def vfit1d(x,data) :
    """ Do a 1D fit """
    fit_p = fitting.LevMarLSQFitter()
    p_init = models.Polynomial1D(degree=1)
    pfit = fit_p(p_init, x, data)
    return pfit
   
def vfit2d(x,y,data) :
    """ Do a 2D fit """
    fit_p = fitting.LevMarLSQFitter()
    p_init = models.Polynomial2D(degree=1)
    pfit = fit_p(p_init, x, y, data)
    return pfit

def read(file) :
    data=fits.open(file)[1].data
    vmicro = data['fparam'][:,2]
    vmacro = data['fparam'][:,7]
    teff = data['fparam'][:,0]
    logg = data['fparam'][:,1]
    mh = data['fparam'][:,3]
    try :
       meanfib = data['meanfib']
    except :
       meanfib = 0.*vmicro
    return teff, logg, mh, vmicro, vmacro, meanfib

def fit_vmacro(file) :
    """ Fit microturbulence relation """
    teff, logg, mh, vmicro, vmacro, meanfib = read(file)
    vrange = [1,15]
    gd = np.where((mh>-1) & (logg < 3.8) & (10.**vmacro < vrange[1]))
    gd=gd[0]
    print 'len(gd)', len(gd)
    pdb.set_trace()
    fit1d = vfit1d(logg[gd], vmacro[gd])
    fit2d = vfit2d(teff[gd], logg[gd], vmacro[gd])
    plot(teff, logg, mh, meanfib, vmacro, vrange, fit1d, fit2d, vt='vmacro')
    return fit1d, fit2d

def fit_vmicro(file) :
    """ Fit microturbulence relation """
    teff, logg, mh, vmicro, vmacro, meanfib = read(file)
    vrange = [0,4]
    gd = np.where((mh>-1) & (logg < 3.8) & (10.**vmicro < vrange[1]))
    fit1d = vfit1d(logg[gd], vmicro[gd])
    fit2d = vfit2d(teff[gd], logg[gd], vmicro[gd])
    plot(teff, logg, mh, meanfib, vmicro, vrange, fit1d, fit2d, vt='vmicro')
    return fit1d, fit2d

def plot(teff, logg, mh, meanfib, v, vr, fit1d, fit2d, 
         xr=[5500,3500], yr=[5,0], vt='vmicro') :
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

    ax=fig.add_subplot(3,2,5)
    plots.plotc(ax, logg, 10.**v, mh, xr=[0,5], yr=vr, zr=[-2.5,0.5],
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')

    ax=fig.add_subplot(3,2,6)
    plots.plotc(ax, logg, 10.**v, meanfib, xr=[0,5], yr=vr, zr=[0,300],
                size=10,colorbar=True,xt='log g', yt=vt, zt='<fiber>')

    plt.show()

def litplot(mass=1, feh=0) : 
#def main(argv) :

#    parser = argparse.ArgumentParser()
#    parser.add_argument("--mass",default=1.)
#    parser.add_argument("--feh",default=0.)
#    args=parser.parse_args()
#    print 'Using mass: ', args.mass
#    mass=float(args.mass)
#    feh=float(args.feh)
#    print 'Using feh: ', feh
#    mass = 0.9

    # setup plots
    mpl.rcParams['xtick.labelsize'] = 8
    mpl.rcParams['ytick.labelsize'] = 8
    p=plt.figure()

    # get isochrone data
    for iso in [0.,-1.,-2.] :
        if iso >=0 : 
            c='p' 
        else: 
            c='m'
        cfeh='z'+c+'{:02d}'.format(int(round(abs(iso)*10)))
        print cfeh
        a=ascii.read(os.getenv('ISOCHRONE_DIR')+'/'+cfeh+'.dat')
        age=a['col2']
        j=np.where((age == 9.6) | (age == 9.0) | (age == 10.0))
        logl=a['col5'][j]
        iso_te=10.**a['col6'][j]
        iso_logg=a['col7'][j]
        iso_feh=iso_logg*0.+iso

        ax=p.add_subplot(221)
        lw=0
        vm=vmacro_massarotti(iso_te,iso_logg,iso_feh)
        plots.plotc(ax,iso_te,iso_logg,vm,xr=[7500,2500],yr=[5,-1.],zr=[0,10],linewidth=lw,size=15)
        ax=p.add_subplot(222)
        vm=vmacro_shetrone(iso_te,iso_logg,iso_feh)
        plots.plotc(ax,iso_te,iso_logg,vm,xr=[7500,2500],yr=[5,-1.],zr=[0,10],linewidth=lw,size=15)
        ax=p.add_subplot(223)
        vm=vmacro_bergemann(iso_te,iso_logg,iso_feh)
        plots.plotc(ax,iso_te,iso_logg,vm,xr=[7500,2500],yr=[5,-1.],zr=[0,10],linewidth=lw,size=15)
        ax=p.add_subplot(224)
        vm=vmacro_bergemann_new(iso_te,iso_logg,iso_feh)
        plots.plotc(ax,iso_te,iso_logg,vm,xr=[7500,2500],yr=[5,-1.],zr=[0,10],linewidth=lw,size=15)


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
    ax=p.add_subplot(221)
    ax.set_title('Massarotti vmacro')
    axim=ax.imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500,2500.,5.,-1.],aspect='auto',origin='upper')
    p.colorbar(axim)

    a=vmacro_shetrone(teff,logl,feh)
    ax=p.add_subplot(222)
    ax.set_title('Shetrone vmacro')
    axim=ax.imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500,2500.,5.,-1.],aspect='auto',origin='upper')
    p.colorbar(axim)
    
    a=vmacro_bergemann(teff,logg,feh)
    ax=p.add_subplot(223)
    ax.set_title('Bergemann vmacro')
    axim=ax.imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500.,2500.,5.,-1.],aspect='auto',origin='upper')
    p.colorbar(axim)

    a=vmacro_bergemann_new(teff,logg,feh)
    ax=p.add_subplot(224)
    ax.set_title('Bergemann vmacro new')
    axim=ax.imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500.,2500.,5.,-1.],aspect='auto',origin='upper')
    p.colorbar(axim)

    p.savefig('vmacro_'+str(mass)+'.jpg')

if __name__ == '__main__' :
    main(sys.argv)

