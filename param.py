import matplotlib.pyplot as plt
import numpy as np
from holtz.apogee import apload
from holtz.apogee import apselect
from holtz.tools import plots
from holtz.tools import html
import pdb
from astropy.io import fits
from astropy.io import ascii

reload(apselect)
#reload(galmodel)
reload(plots)

def read(file=None) :
    '''
    Read allStar file

    Returns:
       allstar structure, elem, elemtoh
    '''
    apload.dr13()
    c=apload.allStar()[3].data
    if file is None :
        a=apload.allStar()[1].data
    else :
        print('reading', file)
        a=fits.open(file)[1].data
    elem=c['ELEM_SYMBOL'][0]
    elemtoh=c['ELEMTOH'][0]
    return a, elem, elemtoh

def plot(a,sn=[0,1000]) :
    '''
    plots a bunch of parameter-related things!

    Args:
        a  : allstar Structure

    Keyword args:
        sn=[snlow,snhigh] : range of S/N to accept (default=[0,1000])

    '''
    gd=apselect.select(a,badval='STAR_BAD',sn=sn,raw=True)
    gdgiant=apselect.select(a,badval='STAR_BAD',sn=sn,raw=True,giants=True)

    # get the indices for different grids, and for stars near solar metallicity
    fgrid=apselect.select(a[gd],grid='F',raw=True)
    gkgrid=apselect.select(a[gd],grid='GK',raw=True)
    mgrid=apselect.select(a[gd],grid='M',raw=True)
    solar=apselect.select(a[gd],mh=[-0.1,0.1],raw=True)
    rgb=apselect.select(a[gd],rgb=True,giants=True,raw=True)
    rc=apselect.select(a[gd],rc=True,giants=True,raw=True)
    inter=apselect.select(a[gd],inter=True,giants=True,raw=True)

    ytit=[]
    files=[]

    name='param'
    print(name)
    fname = 'param/'+name
    # loop over plots
    xtit = []
    file=[]
    for iplot in range(0,5) :
          size = 8
          if iplot == 0 :
              # log g vs Teff
              x = a['PARAM'][gd,0]
              xr = [8000,2500]
              xr = [4000,3500]
              xt= 'Teff (cal)'
              y = a['PARAM'][gd,1]
              yr=[5,-1]
              yr=[1.5,0.0]
              yt = 'log g (cal)'
              z = a['PARAM'][gd,3]
              zr = [-2.5,0.5]
              zt='[M/H]'
              xtit.append('calibrated ')
              size = 12
          elif iplot == 1 :
              # log g vs Teff raw
              x = a['FPARAM'][gd,0]
              xr = [8000,2500]
              xt= 'Teff (raw)'
              y = a['FPARAM'][gd,1]
              yr=[5,-1]
              yt = 'log g (raw)'
              z = a['FPARAM'][gd,3]
              zr = [-2.5,0.5]
              zt='[M/H] (raw)'
              size=1
              xtit.append('raw')
          elif iplot == 2 :
              # chi2 vs Teff
              x = a['PARAM'][gd,0]
              xr = [8000,2500]
              xt= 'Teff'
              y = a['ASPCAP_CHI2'][gd]
              yr=[0,50.]
              yt = 'ASPCAP_CHI2'
              #z = a['ELEM'][gd,ifeh]
              z = a['FPARAM'][gd,3]
              zr = [-2.5,0.5]
              zt='[M/H]'
              xtit.append('Chi^2')
          elif iplot == 3 :
              # vmicro vs log g
              x = a['FPARAM'][gd,1]
              xr = [-1,5]
              xt= 'log g'
              y = 10.**a['PARAM'][gd,2]
              yr=[0.,8]
              yt = 'vmicro'
              z = a['PARAM'][gd,3]
              zr = [-2.5,0.5]
              zt='[M/H]'
              xtit.append('vmicro')
          elif iplot == 4 :
              # vrot/vmacro vs log g
              x = a['FPARAM'][gd,1]
              xr = [-1,5]
              xt= 'log g'
              y = 10.**a['PARAM'][gd,7]
              yr=[0.,8]
              yt = 'vmacro/rotation'
              z = a['PARAM'][gd,3]
              zr = [-2.5,0.5]
              zt='[M/H]'
              xtit.append('vmacro')

          fig=plt.figure(figsize=(10,8))
          ax=fig.add_subplot(111)
          if len(x) > 0 :
                if len(fgrid) > 0 :
                    plots.plotc(ax,x[fgrid],y[fgrid],z[fgrid],xr=xr,yr=yr,zr=zr,colorbar=False,size=size,marker='s',yt=yt,xt=xt,zt=zt)
                if len(gkgrid) > 0 :
                    plots.plotc(ax,x[gkgrid],y[gkgrid],z[gkgrid],xr=xr,yr=yr,zr=zr,size=size,marker='o',yt=yt,xt=xt,zt=zt)
                    #plots.plotc(ax,x[gkgrid],y[gkgrid],'g',xr=xr,yr=yr,size=8,marker='o',yt=yt,xt=xt,zt=zt)
                if len(mgrid) > 0 :
                    plots.plotc(ax,x[mgrid],y[mgrid],z[mgrid],xr=xr,yr=yr,zr=zr,size=size,marker='^',yt=yt,xt=xt,zt=zt)
                    #plots.plotc(ax,x[mgrid],y[mgrid],'r',xr=xr,yr=yr,size=8,marker='^',yt=yt,xt=xt,zt=zt)
          if iplot == 0 or iplot ==1 :
              plots.plotl(ax,[7000.,7000.],[4.,0.])
              plots.plotl(ax,[7000.,4800.],[4.,4.])
              plots.plotl(ax,[4800.,3500.],[4.,2.])
              plots.plotl(ax,[7000.,7000.],[3.75,0.])
              plots.plotl(ax,[7000.,4800.],[3.75,3.75])
              plots.plotl(ax,[4800.,3500.],[3.75,1.75])

          plt.savefig(fname+'{:1d}.jpg'.format(iplot))
          plt.close()
          file.append(name+'{:1d}.jpg'.format(iplot))

    files.append(file)
    ytit.append(name)
    file=[]
    xtit= []
    for iplot in range(5,10) :
          if iplot == 5 :
              # Delta T vs Teff by class
              x = a['FPARAM'][gdgiant,0]
              xr = [3000,4500]
              xt= 'Teff'
              y=[]
              for iclass in [3,4,5,6,11,12,13,14,15,16,17,18] :
                print(iclass)
                y.append(a['FPARAM_CLASS'][gdgiant,iclass,0]-a['FPARAM'][gdgiant,0])
              yr=[-100.,100.]
              yt = 'Delta T'
              z = a['FPARAM'][gdgiant,3]
              zr = [-2.5,0.5]
              zt='[M/H]'
              xtit.append('Delta Teff')
          elif iplot == 6 :
              # Delta logg vs Teff by class
              x = a['FPARAM'][gdgiant,0]
              xr = [3000,4500]
              xt= 'Teff'
              y=[]
              for iclass in [3,4,5,6,11,12,13,14,15,16,17,18] :
                y.append(a['FPARAM_CLASS'][gdgiant,iclass,1]-a['FPARAM'][gdgiant,1])
              yr=[-0.3,0.3]
              yt = 'Delta logg'
              z = a['FPARAM'][gdgiant,3]
              zr = [-2.5,0.5]
              zt='[M/H]'
              xtit.append('Delta logg')
          elif iplot == 7 :
              # Delta [M/H] vs Teff by class
              x = a['FPARAM'][gdgiant,0]
              xr = [3000,4500]
              xt= 'Teff'
              y=[]
              for iclass in [3,4,5,6,11,12,13,14,15,16,17,18] :
                y.append(a['FPARAM_CLASS'][gdgiant,iclass,3]-a['FPARAM'][gdgiant,3])
              yr=[-0.3,0.3]
              yt = 'Delta [M/H]'
              z = a['FPARAM'][gdgiant,3]
              zr = [-2.5,0.5]
              zt='[M/H]'
              xtit.append('Delta [M/H]')
          elif iplot == 8 :
              # Delta [alpha/M] vs Teff by class
              x = a['FPARAM'][gdgiant,0]
              xr = [3000,4500]
              xt= 'Teff'
              y=[]
              for iclass in [3,4,5,6,11,12,13,14,15,16,17,18] :
                y.append(a['FPARAM_CLASS'][gdgiant,iclass,6]-a['FPARAM'][gdgiant,6])
              yr=[-0.5,0.5]
              yt = 'Delta [alpha/H]'
              z = a['FPARAM'][gdgiant,3]
              zr = [-2.5,0.5]
              zt='[M/H]'
              xtit.append('Delta [alpha/H]')
          elif iplot == 9 :
              # Delta logg vs Teff by class
              x = a['FPARAM'][gdgiant,0]
              xr = [3000,4500]
              xt= 'Teff'
              y=[]
              z=[]
              for iclass in [3,4,5,6,11,12,13,14,15,16,17,18] :
                y.append(a['FPARAM_CLASS'][gdgiant,iclass,1]-a['FPARAM'][gdgiant,1])
                z.append(a['FPARAM_CLASS'][gdgiant,iclass,0]-a['FPARAM'][gdgiant,0])
              yr=[-0.3,0.3]
              yt = 'Delta logg'
              zr = [-50,50]
              zt='Delta Teff'
              xtit.append('Delta logg')
          fig=plt.figure(figsize=(10,8))
          ax=fig.add_subplot(111)
          if len(x) > 0 :
              if iplot < 9 :
                for jj in range(0,len(y)) :
                  plots.plotc(ax,x,y[jj],z,xr=xr,yr=yr,zr=zr,colorbar=False,size=10,marker='s',yt=yt,xt=xt,zt=zt)
              else :
                for jj in range(0,len(y)) :
                  plots.plotc(ax,x,y[jj],z[jj],xr=xr,yr=yr,zr=zr,colorbar=False,size=10,marker='s',yt=yt,xt=xt,zt=zt)
          plt.savefig(fname+'{:1d}.jpg'.format(iplot))
          plt.close()
          file.append(name+'{:1d}.jpg'.format(iplot))
    files.append(file)
    ytit.append(name)

    file = []
    xtit= []
    for iplot in range(10,15) :
          x = a['FPARAM'][gd,1]
          xr = [-1,5]
          xt= 'log g'
          y = a['PARAM'][gd,1]-a['FPARAM'][gd,1]
          yr=[-1,1]
          yt = 'delta log g'
          z = a['PARAM'][gd,3]
          zr = [-2.5,0.5]
          zt='[M/H]'
          xtit.append('calibrated ')

          fig=plt.figure(figsize=(10,8))
          ax=fig.add_subplot(111)
          if iplot == 10 :
              plots.plotc(ax,x[rgb],y[rgb],z[rgb],xr=xr,yr=yr,zr=zr,colorbar=False,size=8,marker='s',yt=yt,xt=xt,zt=zt)
          elif iplot == 11 :
              plots.plotc(ax,x[rc],y[rc],z[rc],xr=xr,yr=yr,zr=zr,colorbar=False,size=8,marker='s',yt=yt,xt=xt,zt=zt)
          elif iplot == 12 :
              plots.plotc(ax,x[inter],y[inter],z[inter],xr=xr,yr=yr,zr=zr,colorbar=False,size=8,marker='s',yt=yt,xt=xt,zt=zt)
          plt.savefig(fname+'{:1d}.jpg'.format(iplot))
          plt.close()
          file.append(name+'{:1d}.jpg'.format(iplot))
    files.append(file)
    ytit.append(name)
   
    html.htmltab(files,ytitle=ytit,file='param/'+'param.html',xtitle=xtit)

def alpha(a) :
    '''
    Plot difference between [alpha/M] and individual alpha abundances
    '''
    gd=apselect.select(a,badval='STAR_BAD')

    fig,ax=plots.multi(1,6,hspace=0.001)

    plots.plotc(ax[0],a['TEFF'][gd],a['O_FE'][gd]-a['ALPHA_M'][gd],a['M_H'][gd],xr=[3000,6000],yr=[-0.5,0.5],zr=[-2,0.5],yt='[O/FE]-[alpha/M]')
    plots.plotc(ax[1],a['TEFF'][gd],a['MG_FE'][gd]-a['ALPHA_M'][gd],a['M_H'][gd],xr=[3000,6000],yr=[-0.5,0.5],zr=[-2,0.5],yt='[Mg/FE]-[alpha/M]')
    plots.plotc(ax[2],a['TEFF'][gd],a['SI_FE'][gd]-a['ALPHA_M'][gd],a['M_H'][gd],xr=[3000,6000],yr=[-0.5,0.5],zr=[-2,0.5],yt='[Si/FE]-[alpha/M]')
    plots.plotc(ax[3],a['TEFF'][gd],a['S_FE'][gd]-a['ALPHA_M'][gd],a['M_H'][gd],xr=[3000,6000],yr=[-0.5,0.5],zr=[-2,0.5],yt='[S/FE]-[alpha/M]')
    plots.plotc(ax[4],a['TEFF'][gd],a['CA_FE'][gd]-a['ALPHA_M'][gd],a['M_H'][gd],xr=[3000,6000],yr=[-0.5,0.5],zr=[-2,0.5],yt='[Ca/FE]-[alpha/M]')
    plots.plotc(ax[5],a['TEFF'][gd],a['TI_FE'][gd]-a['ALPHA_M'][gd],a['M_H'][gd],xr=[3000,6000],yr=[-0.5,0.5],zr=[-2,0.5],yt='[Ti/FE]-[alpha/M]',xt='Teff')
    

def main() :
    ''' 
    Make series of plots and web pages for each calibration "type" 
    '''    

    a,e,etoh = read('allStar-l30e.2.fits')
    #a,e,etoh = read('allStar-testcal.fits')
    plot(a)

if __name__ == '__main__' :
    main()
