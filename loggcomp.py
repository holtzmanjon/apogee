from sdss.apogee import apload
from sdss.apogee import apselect
from holtz.tools import match
from holtz.gal import isochrones
import astropy
from astropy.io import fits
from holtz.tools import plots
import pdb
import matplotlib.pyplot as plt
import numpy as np
import os

os.environ['IDLWRAP_DIR']='/home/holtz/sdss4/idlwrap'

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

fig=plt.figure()
ax1=fig.add_subplot(231)
ax2=fig.add_subplot(232)
ax3=fig.add_subplot(233)
ax4=fig.add_subplot(234)
ax5=fig.add_subplot(235)
ax6=fig.add_subplot(236)

clust=apselect.clustdata()
for cluster in ['M92','M15','M53','M2','M13','M3','M5'] :
    i=np.where(clust.name == cluster)
    dist=clust[i].dist*1000.
    mh=clust[i].mh
    mass=0.85
    j=apselect.clustmember(dr12,cluster,raw=True)
    lum=10.**(-0.4*(dr12['H'][j]+isochrones.bc(dr12['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
    logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*dr12['FPARAM'][j,0]**4/lum)
    plots.plotc(ax1,dr12['FPARAM'][j,3]*0+mh,dr12['FPARAM'][j,1]-logg,dr12['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],yt='ASPCAP-physical log g',label=[-2.5,1.5,'DR2 raw'])
    plots.plotc(ax4,dr12['PARAM'][j,3]*0+mh,dr12['PARAM'][j,1]-logg,dr12['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-physical log g',label=[-2.5,1.5,'DR12 cal'])

    j=apselect.clustmember(dr13,cluster,raw=True)
    lum=10.**(-0.4*(dr13['H'][j]+isochrones.bc(dr13['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
    logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*dr13['FPARAM'][j,0]**4/lum)
    plots.plotc(ax2,dr13['FPARAM'][j,3]*0+mh,dr13['FPARAM'][j,1]-logg,dr13['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[-2.5,1.5,'DR13 raw'])
    plots.plotc(ax5,dr13['PARAM'][j,3]*0+mh,dr13['PARAM'][j,1]-logg,dr13['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[-2.5,1.5,'DR13 cal'],xt='[M/H]')

    j=apselect.clustmember(l30i,cluster,raw=True)
    lum=10.**(-0.4*(l30i['H'][j]+isochrones.bc(l30i['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
    logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*l30i['FPARAM'][j,0]**4/lum)
    plots.plotc(ax3,l30i['FPARAM'][j,3]*0+mh,l30i['FPARAM'][j,1]-logg,l30i['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[-2.5,1.5,'l30i raw'])
    plots.plotc(ax6,l30i['PARAM'][j,3]*0+mh,l30i['PARAM'][j,1]-logg,l30i['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[-2.5,1.5,'l30i cal'],xt='[M/H]')
    
    fig.subplots_adjust(hspace=0.001,wspace=0.001)
    ticklabels = ax2.get_yticklabels() + ax3.get_yticklabels() + ax5.get_yticklabels() + ax6.get_yticklabels() + ax1.get_xticklabels() + ax2.get_xticklabels() + ax3.get_xticklabels()
    plt.setp(ticklabels, visible=False)

    plt.show()

#plt.tight_layout()
pdb.set_trace()
ax1.cla()
ax2.cla()
ax3.cla()
ax4.cla()
ax5.cla()
ax6.cla()


i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
plots.plotc(ax1,dr13['FPARAM'][i1,3],dr13['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr13['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[-2.,1.5,'DR13 raw'])
plots.plotc(ax4,dr13['PARAM'][i1,3],dr13['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr13['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-Kepler log g',label=[-2.,1.5,'DR13 cal'])

i1,i2=match.match(dr12['APOGEE_ID'],apokasc['2MASS_ID'])
plots.plotc(ax2,dr12['FPARAM'][i1,3],dr12['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr12['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],label=[-2,1.5,'DR12 raw'])
plots.plotc(ax5,dr12['PARAM'][i1,3],dr12['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr12['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[-2,1.5,'DR12 cal'])

i1,i2=match.match(l30i['APOGEE_ID'],apokasc['2MASS_ID'])
plots.plotc(ax3,l30i['FPARAM'][i1,3],l30i['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],l30i['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],label=[-2,1.5,'l30i raw'])
plots.plotc(ax6,l30i['PARAM'][i1,3],l30i['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],l30i['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[-2,1.5,'l30i cal'])

plt.subplots_adjust(hspace=0.001,wspace=0.001)
ticklabels = ax2.get_yticklabels() + ax3.get_yticklabels() + ax5.get_yticklabels() + ax6.get_yticklabels() + ax1.get_xticklabels() + ax2.get_xticklabels() + ax3.get_xticklabels()
plt.setp(ticklabels, visible=False)
#plt.tight_layout()
plt.show()
