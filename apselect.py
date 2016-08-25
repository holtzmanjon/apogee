import numpy as np
import pdb
import copy
from astropy.io import ascii
import os
from holtz.tools import plots
import matplotlib.pyplot as plt

def select(data,badval=None,logg=[-1,10],teff=[0,10000],mh=[-100.,100.],alpha=[-100.,100.],sn=[0,1000], raw=False, glon=[0,360],glat=[-90,90],grid=None,field=None,giants=None, dwarfs=None,rgb=None, rc=None,inter=None, id=None, redid=None) :
    '''  
    Return indices of requested subsamples from input allStar structure 

    Args:
        data : allStar structure

    Keyword args:
        badval (char) : ASPCAPFLAG value or list of values to reject
        logg  (float[2]) : min, max log g to accept (default [-1,10])
        teff  (float[2]) : min, max Teff to accept (default [0,10000])
        mh    (float[2]) : min, max [M/H] to accept (default [-100,100])
        alpha (float[2]) : min, max [alpha/H] to accept (default [-100,100])
        sn    (float[2]) : min, max S/N to accept (default [0,1000])
        glon  (float[2]) : min, max GLON to accept (default [0,360])
        glat  (float[2]) : min, max GLAT to accept (default [-90,90])
        grid  (char)     : ASPCAP_CLASS to accept, if specified (default None)
        field (char)     : FIELD to accept, if specified (default None)
        giants (bool)    : ASPCAP giants if true, per line in HR diagram, else all (default = None)
        dwarfs (bool)    : ASPCAP dwarfs if true, per line in HR diagram, else all (default = None)
        rgb    (bool)    : ASPCAP RGB per Teff, [M/H], [C/N] criterion (default = None)
        rc     (bool)    : ASPCAP RC per Teff, [M/H], [C/N] criterion (default = None)
        inter  (bool)    : ASPCAP "intermediate" RGB/RC per Teff, [M/H], [C/N] criterion (default = None)
        raw   (bool)     : specifies raw (FPARAM) values rather than calibrated (PARAM) (default=False)
        id    (char)     : APOGEE_ID to return (ignores other constraints)
        redid (char)     : REDUCTION_ID to return (ignores other constraints)
 
    '''

    if id is not None :
        gd = np.where(np.core.defchararray.find(data.APOGEE_ID,id) >=0 )
        return gd[0]
 
    if redid is not None :
        gd = np.where(np.core.defchararray.find(data.REDUCTION_ID,redid) >=0 )
        return gd[0]
    
    if alpha is None :
        alpha=[-100,100]    
    if mh is None :
        mh=[-100,100]    
    if glon is None :
        glon=[0,360]
    if glat is None :
        glat=[-90,90]
 
    if type(badval) is str :
        badval = [badval]
    badbits = 0
    if badval is not None : 
        for val in badval :
            badbits = badbits and flag.aspcapflagval(val)
   
    if raw :
        param='FPARAM'
    else :
        param='PARAM'
    t = data[param][:,0] 
    g = data[param][:,1] 
    m = data[param][:,3] 
    a = data[param][:,6] 

    if giants is not None or dwarfs is not None :
        startype = (g < 2./1300.*(t-3500.)+2) & (g < 4) & (t < 7000)
        if dwarfs is not None :
           startype = np.logical_not(startype)
    else  :
        # take all stars
        startype = g < 100

    # define dt from Bovy ridgeline as f([M/H])
    dt=data['FPARAM'][:,0]-(4468+(data['FPARAM'][:,1]-2.5)/0.0018 - 382.5*data['FPARAM'][:,3])
    cn=data['FPARAM'][:,4]-data['FPARAM'][:,5]
    if rgb :
       startype = np.logical_and(startype,
                  (dt<0) | (data['FPARAM'][:,1] < 2.3) |
                  ((dt>0) & (dt<100) & (cn<-0.1-0.005*dt)))
    if rc :
       startype = np.logical_and(startype,
                  (data['FPARAM'][:,1] > 2.3) &
                  ((dt>100) | ((dt>0) & (cn>-0.075-0.0025*dt)) ))
    if inter :
       startype = np.logical_and(startype,
                  (data['FPARAM'][:,1] > 2.3) & (dt>0) & (dt<100) &
                  (cn>-0.1-0.005*dt) & (cn<-0.075-0.0025*dt))

    try :
       bad = data.ASPCAPFLAG & badbits
    except :
       bad = np.zeros(len(data),dtype=np.int8)

    try :
       snr = data.SNR
    except :
       snr = data['SNR_2']

    gd = np.where((bad == 0)  &
         (startype) &
         (t > teff[0]) & (t < teff[1])  &
         (g > logg[0]) & (g < logg[1])  &
         (m > mh[0]) & (m < mh[1])  &
         (a > alpha[0]) & (a < alpha[1])  &
         (data.GLON > glon[0]) & (data.GLON < glon[1])  &
         (data.GLAT > glat[0]) & (data.GLAT < glat[1])  &
         (snr > sn[0]) & (snr < sn[1])  
         )[0]

    if grid is not None :
        gdclass = np.where(np.core.defchararray.find(data.ASPCAP_CLASS[gd],grid) >=0 )[0]
        gd=gd[gdclass]

    if field is not None :
        gdfield = np.where(np.core.defchararray.strip(data.FIELD[gd]) ==field)[0]
        gd=gd[gdfield]

    return gd

def clustdata() :
    """
    Returns structure containing cluster data
    """

    clust=['M92','M15','M53','N5466','N4147',
        'M2','M13','M3','M5','M12','M107',
        'M71','N2243','Be29', 'N2158','M35','N2420',
        'N188','M67','N7789','Pleiades','N6819',
        'N6791']
    out = np.recarray(len(clust),dtype=[
                       ('name','S24'),
                       ('field','S24'),
                       ('rv','f4'),
                       ('drv','f4'),
                       ('mh','f4'),
                       ('dist','f4'),
                       ('age','f4'),
                       ('giant_mass','f4'),
                       ('ra','f4'),
                       ('dec','f4'),
                       ('rad','f4'),
                       ('ebv','f4'),
                       ])
    out['ebv']=[0.02,0.1,0.02,0.,0.02,
                0.06,0.02,0.01,0.03,0.19,0.33,
                0.25,0.051,0.157,0.36,0.262,0.05,
                0.08,0.04,0.217,0.03,0.16,
                0.122]
    out['name']=clust
    out['field']=['M92','M15','M53','N5466','N4147',
                'M2','M13','M3','M5PAL5','M12','M107',
                'M71','N2243','198+08', 'M35N2158','M35N2158','N2420',
                'N188','M67','N7789','Pleiades','N6819',
                'N6791']
    out['rv']=[ -118.517, -107.508, -61.5988, 106.883, 183., 
              -3.74874, -246.589, -145.525, 54.8727, -41.4000, -35.2638,
              -23.5504, 60., 25., 26.6529, -7.02990, 74.3025, 
              -42.1051, 34.0525, -55.0546, 5.59228, 2.55553, 
              -47.4558]
    out['drv']=[ 10., 10., 10., 10., 10., 
               10., 10., 10., 10., 10., 8., 
               8., 10., 10., 10., 10., 10., 
               6., 6., 6., 10., 6., 
               6.]
    out['mh']=[-2.35,-2.33,-2.06,-2.01,-1.78,
             -1.66,-1.58,-1.50,-1.33,-1.37,-1.03,
             -0.82,-0.35,-0.44,-0.21,-0.14,-0.13,
              0.04, 0.06,0.09,0.03,0.16,
              0.37]
    out['dist']=[8.3, 10.4, 17.9, 16.0, 19.3, 
               12.59, 7.1, 10.2, 7.5, 6.33, 6.4, 
               4.0, 4.45, 14.87, 5.06, 0., 2.44, 
               2.04, 0.907, 2.33, 0., 2.36, 
               4.09]
    out['age']=[12., 12., 12., 12., 12., 
               12., 12., 12., 12., 12., 12., 
               10.0, 4.5, 1.1, 1.0, 0.5, 1.1, 
               4.26, 6., 1.7, 0.5, 1.5, 
               4.4]
    out['giant_mass'] = [ 0.85, 0.85, 0.85, 0.85, 0.85,
                          0.85, 0.85, 0.85, 0.85, 0.85, 0.85,
                          0.85, -9999.,  -9999.,  -9999., -9999., 1.6,
                          -9999., 1.36, -9999., -9999., 1.63,
                          1.15 ]
    out['ra']=[259.27917,322.4929,198.2292,211.3625,182.5262,
             323.3625,250.42083,205.55,229.6375,251.80908,248.1333,
             298.4438,97.3917,103.325,91.8542,91.854,114.5958,
             12.1083,132.825,359.35,56.75,295.325,
             290.22083]
    out['dec']=[43.1358,12.1669,18.1806,28.5344,18.5425,
              -0.82325,36.4611,28.3772,2.0811,-1.94853,-13.05361,
              18.7792,-31.2833,16.9167,24.0967,24.097,21.5733,
              85.255,11.8,56.7083,24.11667,40.1867,
              37.77167]
    out['rad']=[20.,18.,13.,11.,4.,
             16.,20.,24.,23.,16.,13.,
             7.,5.,5.,5.,10.,6.,
             16.,45.,16.,110.,12.,
             16.]

    return out.view(np.recarray)


def clustmember(data,cluster,logg=[-1,3.8],te=[3800,5500],raw=False,firstgen=False,plot=False,hard=None) :

    clust=clustdata()
    ic = np.where( np.core.defchararray.strip(clust.name) == cluster)[0]
    if len(ic) == 0 :
        print 'no cluster found: ',cluster
        return []

    # adjust ra for wraparound if needed
    ra=copy.copy(data['ra'])
    if clust[ic].ra > 300 :
        j=np.where(data['ra'] < 180)[0]
        ra[j]+=360
    if clust[ic].ra < 60 :
        j=np.where(data['ra'] > 180)[0]
        ra[j]-=360

    # select by location relative to cluster
    jc=np.where((np.abs(ra-clust[ic].ra)*np.cos(clust[ic].dec*np.pi/180.) < clust[ic].rad/60.) & 
                (np.abs(data['dec']-clust[ic].dec) < clust[ic].rad/60.))[0]
    if len(jc) > 0 :
        j=np.where( ((ra[jc]-clust[ic].ra)*np.cos(clust[ic].dec*np.pi/180.))**2+ 
                     (data[jc]['dec']-clust[ic].dec)**2 < (clust[ic].rad/60.)**2)[0]
        jc=jc[j]
    else :
        print 'no stars after location criterion'
        jc=[]
    if plot :
        jf=np.where((np.abs(ra-clust[ic].ra)*np.cos(clust[ic].dec*np.pi/180.) < 1.5) & 
                (np.abs(data['dec']-clust[ic].dec) < 1.5))[0]
        fig,ax=plots.multi(1,1)
        plots.plotp(ax,ra[jf],data['dec'][jf],color='k',size=20)
        plots.plotp(ax,ra[jc],data['dec'][jc],color='g',size=20)
        if hard is not None :
            print(hard+'/'+clust[ic].name[0]+'_pos.jpg')
            fig.savefig(hard+'/'+clust[ic].name[0]+'_pos.jpg')
        else :
            pdb.set_trace()

    # RV criterion
    try :
        vhelio = data['vhelio']
    except :
        vhelio = data['vhelio_avg']
    j=np.where(np.abs(vhelio[jc]-clust[ic].rv) < clust[ic].drv)[0]
    if plot :
        ax.cla() 
        ax.hist(vhelio[jf],color='k',bins=np.arange(clust[ic].rv-100,clust[ic].rv+100,1.),histtype='step')
        ax.hist(vhelio[jc],color='g',bins=np.arange(clust[ic].rv-100,clust[ic].rv+100,1.),histtype='step')
        if hard is not None :
            fig.savefig(hard+'/'+clust[ic].name[0]+'_rv.jpg')
        else :
            pdb.set_trace()
    if len(j) > 0 :
        jc=jc[j]
    else :
        print 'no stars after RV criterion'
        jc=[]

    # parameters criteria
    if raw :
        param='FPARAM'
    else :
        param='PARAM'
    try :
        j = np.where((data[param][jc,1] >= logg[0]) & (data[param][jc,1] <= logg[1]) &
                     (data[param][jc,0] >= te[0]) & (data[param][jc,0] <= te[1]) )[0]
        if len(j) > 0 :
            jc=jc[j]
        else :
            jc=[]
            print 'no stars after parameters criterion'
    except: pass

    # Remove badstars
    if plot :
        ax.cla()
        plots.plotp(ax,data['J'][jf]-data['K'][jf],data['K'][jf],color='k',size=20,xr=[-0.5,1.5],yr=[15,6],facecolors='none',linewidth=1)
        plots.plotp(ax,data['J'][jc]-data['K'][jc],data['K'][jc],color='g',size=30,xr=[-0.5,1.5],yr=[15,6])
    badstars = open(os.environ['IDLWRAP_DIR']+'/data/badcal.dat')
    bad = []
    for line in badstars :
       bad.append(line.split()[0])
    jc = [x for x in jc if data[x]['APOGEE_ID'] not in bad]

    # remove non firstgen GC stars if requested
    if firstgen :
        gcstars = ascii.read(os.environ['IDLWRAP_DIR']+'/data/gc_szabolcs.dat')
        #bd=np.where(gcstars['pop'] != 1)[0]
        #jc = [x for x in jc if data[x]['APOGEE_ID'] not in gcstars['id'][bd]]
        gd=np.where(gcstars['pop'] == 1)[0]
        jc = [x for x in jc if data[x]['APOGEE_ID'] in gcstars['id'][gd]]

    if plot :
        plots.plotp(ax,data['J'][jc]-data['K'][jc],data['K'][jc],color='b',size=30)
        if hard is not None :
            fig.savefig(hard+'/'+clust[ic].name[0]+'_cmd.jpg')
        else :
            pdb.set_trace()

    return jc
 
