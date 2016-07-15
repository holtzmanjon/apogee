from sdss.apogee import apload
from holtz.tools import match
from astropy.io import fits
from holtz.tools import plots
from holtz.apogee import flag
import pdb
import matplotlib.pyplot as plt
import numpy as np

apokasc=fits.open('APOKASC_cat_v3.6.0.fits')[1].data
logg='LOGG_SYD_SCALING'
#apokasc=fits.open('APOKASC_cat_v3.3.5.fits')[1].data
#logg='LOGGRG'

#j=np.where((apokasc['TEFF_FIT'] < 4000) & (apokasc[logg] > -500))[0]
#j=np.where((apokasc['CONS_EVSTATES'] == 'RGB'))[0]
#pdb.set_trace()
#apokasc=apokasc[j]

apload.dr13()
dr13=apload.allStar()[1].data
apload.aspcap='l30i'
apload.results='l30i'
l30i=apload.allStar()[1].data

fig=plt.figure()
ax1=fig.add_subplot(211)
ax2=fig.add_subplot(212)
#ax3=fig.add_subplot(223)
#ax4=fig.add_subplot(224)

i1,i2=match.match(l30i['APOGEE_ID'],apokasc['2MASS_ID'])
warn=np.where(l30i['ASPCAPFLAG'][i1] & flag.aspcapflagval('ATMOS_HOLE_WARN'))[0]
bad=np.where(l30i['ASPCAPFLAG'][i1] & flag.aspcapflagval('ATMOS_HOLE_BAD'))[0]
#plt.plot(apokasc[logg][i2],l30i['FPARAM'][i1,1],'ro')
#plots.plotc(ax1,l30i['FPARAM'][i1,0],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,3],zr=[-1,0.5],xr=[3500,5000],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'l30i (MARCS)'],colorbar=True,zt='[M/H]')
plots.plotc(ax1,l30i['FPARAM'][i1,0],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,6],zr=[-0.1,0.4],xr=[3500,5000],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'l30i (MARCS)'],colorbar=True,zt='[alpha/M]')
plots.plotp(ax1,l30i['FPARAM'][i1[warn],0],l30i['FPARAM'][i1[warn],1]-apokasc[logg][i2[warn]],color='y')
plots.plotp(ax1,l30i['FPARAM'][i1[bad],0],l30i['FPARAM'][i1[bad],1]-apokasc[logg][i2[bad]],color='r')
#plots.plotc(ax2,l30i['FPARAM'][i1,1],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,0],zr=[3500,5000],xr=[0,5],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[4,0.7,'l30i'],colorbar=True,zt='Teff')

i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
#plt.plot(apokasc[logg][i2],dr13['FPARAM'][i1,1],'bo')
#plt.xlim(1.60,1.15)
#plt.ylim(1.15,1.85)
#plots.plotc(ax2,dr13['FPARAM'][i1,0],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,3],zr=[-1,0.5],xr=[3500,5000],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'dr13 (Kurucz)'],colorbar=True,zt='[M/H]')
plots.plotc(ax2,dr13['FPARAM'][i1,0],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,6],zr=[-0.1,0.4],xr=[3500,5000],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[3600,0.7,'dr13 (Kurucz)'],colorbar=True,zt='[alpha/M]')
##plots.plotc(ax4,dr13['FPARAM'][i1,1],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,0],zr=[3500,5000],xr=[0,5],yr=[-1,1],xt='Kepler log g',yt='ASPCAP-Kepler log g',label=[4,0.7,'dr13'],colorbar=True,zt='Teff')

plt.tight_layout()
plt.show()
