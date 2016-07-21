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

def solar(el) :

    elem=np.array(['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Cr','Mn','Fe','Co','Ni'])
    abund=np.array([8.39,7.78,8.66,6.17,7.53,6.37,7.51,7.14,5.08,6.31,4.90,4.00,5.64,5.39,7.45,4.92,6.23])
    j=np.where(elem == el.strip())[0]
    if len(j) > 0 : 
       return abund[j] 
    else :
       return 0.

def main(mh=False,alpham=False) :
    '''
    
    '''

    apload.dr13()
    a=fits.open('APOGEE_GES_overlap_calib_DR13.fits')[1].data
    all=apselect.select(a,rc=False,rgb=False,badval='STAR_BAD')
    rc=apselect.select(a,rc=True)
    rgb=apselect.select(a,rgb=True)
    unk=np.delete(all,np.append(rc,rgb))

    types=['fts','solar neighborhood','solar']
    size=[40,10,40]
    marker=['o','^','s']
    cbbox=[0.3, 0.95, 0.5, 0.03] 

    if mh or alpham :
      if mh :
        x='FE_H'
        xr=[-2.5,0.75]
      else :
        x='ALPHA_M'
        xr=[-0.5,0.5]
      z='TEFF_2'
      zr=[3500,6000]
    else :
      x='TEFF_2'
      xr=[3500,6000]
      z='M_H'
      zr=[-2,0.5]

    # Teff
    n=3
    name='param'
    fig,ax=plots.multi(1,n,figsize=(8,2+n),hspace=0.001)
    for i in range(n) :
      if i == 0 :
        y1='TEFF_2'
        y2='TEFF_1'
        yr=[-490,490]
      elif i== 1 :
        y1='LOGG_2'
        y2='LOGG_1'
        yr=[-1.5,1.5]
      else :
        y1='ALPHA_M'
        y2='ALPHA_FE'
        yr=[-0.45,0.45]
      map=plots.plotc(ax[i],a[x][unk],a[y1][unk]-a[y2][unk],z=a[z][unk],xt=x,yt=y1+'-'+y2,xr=xr,yr=yr,zr=zr,size=40)
      plots.plotc(ax[i],a[x][rc],a[y1][rc]-a[y2][rc],z=a[z][rc],size=40,marker='^')
      plots.plotc(ax[i],a[x][rgb],a[y1][rgb]-a[y2][rgb],z=a[z][rgb],size=40,marker='s')
      if i == 0 :
        gd=np.where(np.abs(a[y1][all]-a[y2][all])<500.)
      elif i == 1 :
        gd=np.where(np.abs(a[y1][all]-a[y2][all])<2.)
      elif i == 2 :
        gd=np.where(np.abs(a[y1][all]-a[y2][all])<1.)
      ax[i].text(0.05,0.7,'{:<8.2f}'.format(np.nanmean(a[y1][all[gd]]-a[y2][all[gd]])),transform=ax[i].transAxes)
      ax[i].text(0.05,0.6,'{:<8.2f}'.format(np.nanstd(a[y1][all[gd]]-a[y2][all[gd]])),transform=ax[i].transAxes)
    cbaxes = fig.add_axes(cbbox)
    cb = fig.colorbar(map, cax = cbaxes,orientation='horizontal')  
    fig.savefig(name+'.jpg')
    xdata=a[x][all]
    ydata=a[y1][all]-a[y2][all]
    plots._data_x = xdata[np.isfinite(xdata)]
    plots._data_y = ydata[np.isfinite(ydata)]
    plots.event(fig)
    pdb.set_trace()

    # CNO
    n=3
    name='cno'
    fig,ax=plots.multi(1,n,figsize=(8,2+n),hspace=0.001)
    for i in range(n) :
      yr=[-0.45,0.45]
      if i == 0 :
        y1='C_FE'
        y2='C1'
        offset=solar('C')
      elif i== 1 :
        y1='N_FE'
        y2='N2'
        offset=solar('N')
      elif i== 2 :
        y1='O_FE'
        y2='O1'
        offset=solar('O')
      plots.plotc(ax[i],a[x][unk],a[y1][unk]+a['FE_H'][unk]-a[y2][unk]+offset,z=a[z][unk],xt=x,yt=y1+'-'+y2,xr=xr,yr=yr,zr=zr,size=40)
      plots.plotc(ax[i],a[x][rc],a[y1][rc]+a['FE_H'][rc]-a[y2][rc]+offset,z=a[z][rc],size=40,marker='^')
      plots.plotc(ax[i],a[x][rgb],a[y1][rgb]+a['FE_H'][rgb]-a[y2][rgb]+offset,z=a[z][rgb],size=40,marker='s')
      gd=np.where(np.abs(a[y1][all]+a['FE_H'][all]-a[y2][all]+offset)<1.)
      ax[i].text(3600.,0.3,'{:<8.2f}'.format(np.nanmean(a[y1][gd]+a['FE_H'][gd]-a[y2][gd]+offset)))
      ax[i].text(3600.,0.15,'{:<8.2f}'.format(np.nanstd(a[y1][gd]+a['FE_H'][gd]-a[y2][gd]+offset)))
    cbaxes = fig.add_axes(cbbox)
    cb = fig.colorbar(map, cax = cbaxes,orientation='horizontal')  
    fig.savefig(name+'.jpg')
    pdb.set_trace()

    # alpha
    n=5
    name='alpha'
    fig,ax=plots.multi(1,n,figsize=(8,2+n),hspace=0.001)
    for i in range(n) :
      yr=[-0.45,0.45]
      if i == 0 :
        y1='MG_FE'
        y2='MG1'
        offset=solar('Mg')
      elif i== 1 :
        y1='SI_FE'
        y2='SI1'
        offset=solar('Si')
      elif i== 2 :
        y1='S_FE'
        y2='S1'
        offset=solar('S')
      elif i== 3 :
        y1='CA_FE'
        y2='CA1'
        offset=solar('Ca')
      elif i== 4 :
        y1='TI_FE'
        y2='TI1'
        offset=solar('Ti')
      plots.plotc(ax[i],a[x][unk],a[y1][unk]+a['FE_H'][unk]-a[y2][unk]+offset,z=a[z][unk],xt=x,yt=y1+'-'+y2,xr=xr,yr=yr,zr=zr,size=40)
      plots.plotc(ax[i],a[x][rc],a[y1][rc]+a['FE_H'][rc]-a[y2][rc]+offset,z=a[z][rc],size=40,marker='^')
      plots.plotc(ax[i],a[x][rgb],a[y1][rgb]+a['FE_H'][rgb]-a[y2][rgb]+offset,z=a[z][rgb],size=40,marker='s')
      gd=np.where(np.abs(a[y1][all]+a['FE_H'][all]-a[y2][all]+offset)<1.)
      ax[i].text(3600.,0.3,'{:<8.2f}'.format(np.nanmean(a[y1][all[gd]]+a['FE_H'][all[gd]]-a[y2][all[gd]]+offset)))
      ax[i].text(3600.,0.15,'{:<8.2f}'.format(np.nanstd(a[y1][all[gd]]+a['FE_H'][all[gd]]-a[y2][all[gd]]+offset)))
    cbaxes = fig.add_axes(cbbox)
    cb = fig.colorbar(map, cax = cbaxes,orientation='horizontal')  
    fig.savefig(name+'.jpg')
    pdb.set_trace()

    # others
    n=7
    name='nonalpha'
    fig,ax=plots.multi(1,n,figsize=(8,2+n),hspace=0.001)
    el=['Na','Al','K','V','Mn','Fe','Ni']
    ges=['NA1','AL1','K','V1','MN1','FE1','NI1']
    for i in range(n) :
      yr=[-0.45,0.45]
      y1=el[i].upper()+'_FE'
      y2=ges[i]
      offset=solar(el[i])
      if i == 5 :
        plots.plotc(ax[i],a[x][unk],a['FE_H'][unk]-a[y2][unk]+offset,z=a[z][unk],xt=x,yt='FE_H'+'-'+y2,xr=xr,yr=yr,zr=zr,size=40)
        plots.plotc(ax[i],a[x][rc],a['FE_H'][rc]-a[y2][rc]+offset,z=a[z][rc],size=40,marker='^')
        plots.plotc(ax[i],a[x][rgb],a['FE_H'][rgb]-a[y2][rgb]+offset,z=a[z][rgb],size=40,marker='s')
        gd=np.where(np.abs(a['FE_H'][all]-a[y2][all]+offset)<1.)
        ax[i].text(3600.,0.5,'{:<8.2f}'.format(np.nanmean(a['FE_H'][all[gd]]-a[y2][all[gd]]+offset)))
        ax[i].text(3600.,0.25,'{:<8.2f}'.format(np.nanstd(a['FE_H'][all[gd]]-a[y2][all[gd]]+offset)))
      else :
        plots.plotc(ax[i],a[x][unk],a[y1][unk]+a['FE_H'][unk]-a[y2][unk]+offset,z=a[z][unk],xt=x,yt=y1+'-'+y2,xr=xr,yr=yr,zr=zr,size=40)
        plots.plotc(ax[i],a[x][rc],a[y1][rc]+a['FE_H'][rc]-a[y2][rc]+offset,z=a[z][rc],size=40,marker='^')
        plots.plotc(ax[i],a[x][rgb],a[y1][rgb]+a['FE_H'][rgb]-a[y2][rgb]+offset,z=a[z][rgb],size=40,marker='s')
        gd=np.where(np.abs(a[y1][all]+a['FE_H'][all]-a[y2][all]+offset)<1.)
        ax[i].text(3600.,0.3,'{:<8.2f}'.format(np.nanmean(a[y1][all[gd]]+a['FE_H'][all[gd]]-a[y2][all[gd]]+offset)))
        ax[i].text(3600.,0.15,'{:<8.2f}'.format(np.nanstd(a[y1][all[gd]]+a['FE_H'][all[gd]]-a[y2][all[gd]]+offset)))
    cbaxes = fig.add_axes(cbbox)
    cb = fig.colorbar(map, cax = cbaxes,orientation='horizontal')  
    fig.savefig(name+'.jpg')
    pdb.set_trace()




    
