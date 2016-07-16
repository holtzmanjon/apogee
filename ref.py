# routines for comparing ASPCAP results with reference parameters and abundances

import os
import numpy as np
import matplotlib.pyplot as plt
from sdss.apogee import apload
from holtz.tools import match
from holtz.tools import plots
from astropy.io import fits
import pdb

def get_ref(all) :
    """
    Given allStar structure, match by ID to reference catalog
    Args :
           all : allStar structure 
    Returns:
           calindex :  list of indices in cal that have Kepler matches
           ref :       Kepler structure
           refindex :  list of indices in kep that have APOGEE matches
 
    """
    ref=fits.open(os.getenv('IDLWRAP_DIR')+'/data/ref.fits')[1].data

    m1,m2=match.match(all['REDUCTION_ID'].strip(),ref['ID'].strip())
    return m1, ref, m2

def plot(ax,a,r,y,xaxis='TEFF',zaxis='M_H',types=None,size=[40],xr=[3500,6000],zr=[-2,0.5],marker=['o'],yr=None,yt='TEFF') :
    '''
    Plots a series of points with different styles
    '''
    for i in range(len(types)) :
        gd = np.where(r['class'] == types[i])[0]
        print types[i], len(gd), a[xaxis][gd],y[gd]
        sz= size[i] if (len(size) > 1)  else size[0]
        mark=marker[i] if (len(marker) > 1) else marker[0]
        plots.plotc(ax,a[xaxis][gd],y[gd],a[zaxis][gd],xt=xaxis,yt='$\Delta$ '+yt,xr=xr,yr=yr,zr=zr,size=sz,marker=mark)
        ax.text


def main():
    apload.dr13()
    a=apload.allStar()[1].data
    m1,ref,m2 = get_ref(a)
    a=a[m1]
    ref=ref[m2]
    types=['fts','solar neighborhood','solar']
    size=[40,10,40]
    marker=['o','^','s']

    # Teff
    n=2
    fig,ax=plots.multi(1,n,figsize=(8,2+n),hspace=0.001)
    #fig,ax=plt.subplots(n,figsize=(8,2+n))
    #fig.subplots_adjust(hspace=0.001,wspace=0.001)
    #plot(ax[0],a,ref,a['TEFF']-ref['TEFF'],types=types,size=size,marker=marker,yr=[-490,490],yt='Teff')
    #plot(ax[1],a,ref,a['LOGG']-ref['LOGG'],types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='log g')

    plots.plotp(ax[0],a['TEFF'],a['TEFF']-ref['TEFF'],z=a['M_H'],xt='Teff',yt='$\Delta T_{eff}$',xr=[3500,6000],yr=[-490,490],zr=[-2,0.5],
                typeref=ref['class'],types=types,size=size,marker=marker)
    plots.plotp(ax[1],a['TEFF'],a['LOGG']-ref['LOGG'],z=a['M_H'],xt='Teff',yt='$\Delta log g$',xr=[3500,6000],yr=[-0.9,0.9],zr=[-2,0.5],
                typeref=ref['class'],types=types,size=size,marker=marker)
    #ticklabels = ax[0].get_xticklabels() 
    #plt.setp(ticklabels, visible=False)
    pdb.set_trace()

    # CNO
    n=3
    fig,ax=plots.multi(1,n,figsize=(8,2+n),hspace=0.001)
    plots.plotp(ax[0],a['TEFF'],a['C_FE']-(ref['C_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [C/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[1],a['TEFF'],a['N_FE']-(ref['N_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [N/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[2],a['TEFF'],a['O_FE']-(ref['O_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [O/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    #fig.subplots_adjust(hspace=0.001,wspace=0.001)
    #plot(ax[0],a,ref,a['C_FE']-(ref['C_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[C/Fe]')
    #plot(ax[1],a,ref,a['N_FE']-(ref['N_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[N/Fe]')
    #plot(ax[2],a,ref,a['O_FE']-(ref['O_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[O/Fe]')
    #ticklabels = ax[0].get_xticklabels()
    #for i in range(1,n-1) : ticklabels = ticklabels + ax[i].get_xticklabels() 
    #plt.setp(ticklabels, visible=False)
    pdb.set_trace()

    # alpha
    n=5
    fig,ax=plots.multi(1,n,figsize=(8,2+n),hspace=0.001)
    plots.plotp(ax[0],a['TEFF'],a['MG_FE']-(ref['MG_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [Mg/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[1],a['TEFF'],a['SI_FE']-(ref['SI_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [Si/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[2],a['TEFF'],a['S_FE']-(ref['S_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [S/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[3],a['TEFF'],a['CA_FE']-(ref['CA_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [Ca/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[4],a['TEFF'],a['TI_FE']-(ref['TI_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [Ti/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    #fig,ax=plt.subplots(n,figsize=(8,2+n))
    #fig.subplots_adjust(hspace=0.001,wspace=0.001)
    #plot(ax[0],a,ref,a['MG_FE']-(ref['MG_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Mg/Fe]')
    #plot(ax[1],a,ref,a['SI_FE']-(ref['SI_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Si/Fe]')
    #plot(ax[2],a,ref,a['S_FE']-(ref['S_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[S/Fe]')
    #plot(ax[3],a,ref,a['CA_FE']-(ref['CA_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Ca/Fe]')
    #plot(ax[4],a,ref,a['TI_FE']-(ref['TI_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Ti/Fe]')
    #ticklabels = ax[0].get_xticklabels()
    #for i in range(1,n-1) : ticklabels = ticklabels + ax[i].get_xticklabels() 
    #plt.setp(ticklabels, visible=False)
    pdb.set_trace()

    # others
    n=7
    fig,ax=plots.multi(1,n,figsize=(8,2+n),hspace=0.001)
    plots.plotp(ax[0],a['TEFF'],a['NA_FE']-(ref['NA_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [Na/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[1],a['TEFF'],a['AL_FE']-(ref['AL_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [AL/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[2],a['TEFF'],a['K_FE']-(ref['K_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [K/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[3],a['TEFF'],a['V_FE']-(ref['V_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [V/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[4],a['TEFF'],a['MN_FE']-(ref['MN_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [Mn/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[5],a['TEFF'],a['FE_H']-ref['FE_H'],z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [Fe/H]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    plots.plotp(ax[6],a['TEFF'],a['NI_FE']-(ref['NI_H']-ref['FE_H']),z=a['M_H'],typeref=ref['class'],
                xt='Teff',yt='$\Delta [Ni/Fe]$',xr=[3500,6000],yr=[-0.9,0.9],types=types,size=size,marker=marker)
    #fig,ax=plt.subplots(n,figsize=(8,2+n))
    #fig.subplots_adjust(hspace=0.001,wspace=0.001)
    #plot(ax[0],a,ref,a['NA_FE']-(ref['NA_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Na/Fe]')
    #plot(ax[1],a,ref,a['AL_FE']-(ref['AL_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Al/Fe]')
    #plot(ax[2],a,ref,a['K_FE']-(ref['K_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[K/Fe]')
    #plot(ax[3],a,ref,a['V_FE']-(ref['V_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[V/Fe]')
    #plot(ax[4],a,ref,a['MN_FE']-(ref['MN_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Mn/Fe]')
    #plot(ax[5],a,ref,a['FE_H']-ref['FE_H'],types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Fe/H]')
    #plot(ax[6],a,ref,a['NI_FE']-(ref['NI_H']-ref['FE_H']),types=types,size=size,marker=marker,yr=[-0.9,0.9],yt='[Ni/Fe]')
    #ticklabels = ax[0].get_xticklabels()
    #for i in range(1,n-1) : ticklabels = ticklabels + ax[i].get_xticklabels() 
    #plt.setp(ticklabels, visible=False)
    pdb.set_trace()




    
