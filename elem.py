# routines related to individual element calibration for APOGEE/ASPCAP

import matplotlib.pyplot as plt
import numpy as np
from sdss.apogee import apload
from sdss.apogee import apselect
from holtz.tools import plots
from holtz.tools import html
from holtz.tools import fit
import pdb
from astropy.io import fits
from astropy.io import ascii
import esutil
import copy

reload(apselect)
#reload(galmodel)
reload(plots)

def read(file='allStar-testcal.fits') :
    '''
    Read allStar file
    '''
    apload.dr13()
    #a=apload.allStar()[1].data
    c=apload.allStar()[3].data
    #a=fits.open('../dist/allStar+.fits')[1].data
    #x,y,z,r = galmodel.lbd2xyz(a['GLON'],a['GLAT'],a['DISO'][:,2]/1000.)
    #zone=np.where((r>9) & (r<11) & (dt<40))[0]

    a=fits.open(file)[1].data
    #c=fits.open(file)[3].data
    elem=c['ELEM_SYMBOL'][0]
    elemtoh=c['ELEMTOH'][0]

    return a, elem, elemtoh

def arctabun(el) :
    '''
    Define Arcturus abundances, and return requested abundance
    '''
    abun = { "C" : 0.090000, "CI" : 0.09, "N" : 0.400000, "O" : 0.480000, "Na" : 0.210000, "Mg" : 0.370000, "Al" : 0.400000, "Si" : 0.330000, "P" : 0.070000, "S" : 0.350000, "K" : 0.200000, "Ca" : 0.090000, "Sc" : 0.070000, "Ti" : 0.250000, "TiII" : 0.25, "V" : 0.160000, "Cr" : -0.050000, "Mn" : -0.120000, "Fe" : -0.000000, "Co" : 0.040000, "Ni" : 0.030000, "Cu" : -0.050000, "Ge" : 0.000000, "Rb" : 0.000000, "Y" : 0.000000, "Ce" : -0.190000, "Nd" : 0.130000, "M" : 0., "alpha" : 0.3}
    return(abun[el]) 

def optabun(el) :
    '''
    ??? define abundance offsets from some optical analysis ???
    '''
    abun = {"Na" : -0.15, "Mg" :  0.06, "Al" :  0.04, "Si" : -0.21, "Ca" :  0.11, "Ti" : -0.14, "TiII" : 0.08, "V" : -0.15, "Cr" : -0.04, "Mn" : -0.36, "Fe" : 0.06, "Co" : -0.26}
    try :
        return(abun[el]) 
    except :
        return(-9999.)


def refabun(el,dwarf=False) :
   '''
   Return reference abundance: 0 if giant,  Arcturus if not?
   '''
   if dwarf :
       return 0.
   else :
       return arctabun(el)

def plot(a,elem,etoh,dwarf=False,suffix='',gcal=None,dcal=None,glon=None,glat=None,res=None,usemh=False,sn=[200,1000]) :
    '''
    Make a bunch of plots
    '''
    # selection
    #dt=a['FPARAM'][:,0]-(4468+(a['FPARAM'][:,1]-2.5)/0.0018 - 382.5*a['FPARAM'][:,3])
    #gd=apselect.select(a[zone],badval='STAR_BAD',logg=[-1,3.5],sn=[200,1000],teff=[4000,4800])
    #gd=zone[gd]
    if dwarf :
        tit = 'Dwarfs, S/N>200'
        prefix = 'd'+suffix
        tmax=6500
        gd=apselect.select(a,badval='STAR_BAD',sn=sn,raw=True,glon=glon,glat=glat,dwarfs=True)
        etoh[0]=1
        etoh[1]=1
        etoh[2]=1
        ref=apselect.select(a,redid='VESTA')
    else :
        tit = 'Giants, S/N>200'
        prefix = 'g'+suffix
        tmax=6500
        gd=apselect.select(a,badval='STAR_BAD',sn=sn,raw=True,glon=glon,glat=glat,giants=True)
        ref=apselect.select(a,redid='alpha_Boo')
    out = open('elem/'+prefix+'.dat','w')


    # get the indices for different grids, and for stars near solar metallicity
    fgrid=apselect.select(a[gd],grid='F',raw=True)
    gkgrid=apselect.select(a[gd],grid='GK',raw=True)
    mgrid=apselect.select(a[gd],grid='M',raw=True)
    solar=apselect.select(a[gd],mh=[-0.1,0.1],raw=True)

    ifeh=17
    ytit=[]
    files=[]
    # loop over elements
    nelem=len(elem)
    for ielem in range(nelem+2) :
        file=[]
        if ielem < nelem :
            el = elem[ielem].strip()
            #eelem = a['ELEM'][gd,ielem]
            eelem = a['X_M'][gd,ielem]
            felem = a['FELEM'][gd,ielem]
            eelem_err = a['X_M_ERR'][gd,ielem]
            felem_err = a['FELEM_ERR'][gd,ielem]
            tmp=etoh[ielem]
            if ielem > 2 :
                if usemh and etoh[ielem] :
                    #eelem -= a['FPARAM'][gd,3]
                    felem -= a['FPARAM'][gd,3]
                elif not usemh and not etoh[ielem] :
                    #eelem += a['FPARAM'][gd,3]
                    felem += a['FPARAM'][gd,3]
            else :
                giants = apselect.select(a[gd],grid='g_',raw=True)
                if not usemh :
                    #eelem[giants] += a['FPARAM'][gd[giants],3] 
                    felem[giants] += a['FPARAM'][gd[giants],3] 
                dwarfs = apselect.select(a[gd],grid='d_',raw=True)
                if usemh :
                    #eelem[dwarfs] -= a['FPARAM'][gd[dwarfs],3] 
                    felem[dwarfs] -= a['FPARAM'][gd[dwarfs],3] 
        elif ielem == nelem :
            el = 'M'
            eelem = a['PARAM'][gd,0]
            felem = a['FPARAM'][gd,3]
            eelem_err = np.sqrt(a['PARAM_COV'][gd,3,3])
            felem_err = np.sqrt(a['FPARAM_COV'][gd,3,3])
            tmp = 1
        else :
            el = 'alpha'
            eelem = a['PARAM'][gd,6]
            felem = a['FPARAM'][gd,6]
            eelem_err = np.sqrt(a['PARAM_COV'][gd,6,6])
            felem_err = np.sqrt(a['FPARAM_COV'][gd,6,6])
            tmp = 0
            if not usemh :
                eelem += a['FPARAM'][gd,3]
                felem += a['FPARAM'][gd,3]

        if (tmp == 1 and not usemh) or (tmp == 0 and usemh ):
            refoffset=0
        else :
            refoffset=a['FPARAM'][ref,3]
            if usemh :
                refoffset *= -1

        name=prefix+el
        print name
        fname = 'elem/'+name
        # loop over plots
        xtit = []
        for iplot in range(0,8) :
        #for iplot in range(2,3) :
            if iplot == 0 :
              #x = a['ELEM'][gd,ifeh]
              x = a['X_H'][gd,ifeh]
              xr = [-1.5,1.]
              xt= '[Fe/H] (cal)'
              y = eelem
              if not usemh: y-=a['PARAM'][gd,3]
              yr=[-0.25,0.5]
              yt = '['+name+'/M](cal)'
              z = a['FPARAM'][gd,0]
              zr = [3000,tmax]
              zt='Teff'
              xtit.append('calibrated vs [Fe/H]')
            elif iplot == 1 :
              x = a['FELEM'][gd,ifeh]
              xr = [-1.5,1.]
              xt= '[Fe/H] (raw)'
              y = felem
              if not usemh: y-=a['PARAM'][gd,3]
              yr=[-0.25,0.5]
              yt = '['+name+'/M](raw)'
              z = a['FPARAM'][gd,0]
              zr = [3000,tmax]
              zt='Teff'
              xtit.append('raw vs [Fe/H]')
            elif iplot == 2 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt= 'Teff'
              y = eelem
              if not usemh: y-=a['PARAM'][gd,3]
              yr=[-0.25,0.5]
              yt = '['+name+'/M](cal)'
              #z = a['ELEM'][gd,ifeh]
              z = a['X_H'][gd,ifeh]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('calibrated vs Teff')
            elif iplot == 3 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt= 'Teff'
              y = felem
              if not usemh: y-=a['PARAM'][gd,3]
              yr=[-0.25,0.5]
              yt = '['+name+'/M](raw)'
              z = a['FELEM'][gd,ifeh]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('raw vs Teff')
            elif iplot == 4 :
              x = a['FPARAM'][gd,0]
              xr = [3000,tmax]
              xt = 'Teff'
              y = eelem-felem
              yr = [-0.3,0.3]
              yt = 'cal - raw'
              z = a['FELEM'][gd,ifeh]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('calibration')
            elif iplot == 5 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt = 'Teff'
              y = eelem_err
              yr= [0,0.3]
              yt = 'Empirical uncertainty'
              z = a['FELEM'][gd,ifeh]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('empirical uncertainty')
            elif iplot == 6 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt = 'Teff'
              y = felem_err
              yr= [0,0.3]
              yt = 'FERRE uncertainty'
              z = a['FELEM'][gd,ifeh]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('FERRE uncertainty')
            elif iplot == 7 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt = 'Teff'
              if ielem < nelem :
                y = a['ELEM_CHI2'][gd,ielem]
              else :
                y = x*0.
              yr= [0,50]
              yt = 'ELEM_CHI2'
              z = a['FELEM'][gd,ifeh]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('CHI2 from element fit')
    
            fig=plt.figure(figsize=(10,8))
            ax=fig.add_subplot(111)
            if len(x) > 0 :
                if len(fgrid) > 0 :
                    plots.plotc(ax,x[fgrid],y[fgrid],z[fgrid],xr=xr,yr=yr,zr=zr,colorbar=False,size=10,marker='s',yt=yt,xt=xt,zt=zt)
                if len(gkgrid) > 0 :
                    plots.plotc(ax,x[gkgrid],y[gkgrid],z[gkgrid],xr=xr,yr=yr,zr=zr,size=10,marker='o',yt=yt,xt=xt,zt=zt)
                if len(mgrid) > 0 :
                    plots.plotc(ax,x[mgrid],y[mgrid],z[mgrid],xr=xr,yr=yr,zr=zr,size=7,marker='^',yt=yt,xt=xt,zt=zt)
            if (iplot == 0 or iplot ==  2) : 
                if res is not None :
                  clust, = np.where(res['col3'] == ielem)
                  plots.plotp(ax,res['col4'][clust],res['col9'][clust],xr=xr,yr=yr,size=50,marker='o',facecolors='none',linewidth=1)

                # plot the reference star abundance (Arcturus or Vesta)
                if ielem < nelem-2 :
                  #refval = a['ELEM'][ref,ielem]+refoffset
                  #referr = a['ELEM_ERR'][ref,ielem]
                  refval = a['X_M'][ref,ielem]
                  referr = a['X_M_ERR'][ref,ielem]
                elif ielem == nelem-2 :
                  refval = a['PARAM'][ref,3]+refoffset
                  referr = np.sqrt(a['PARAM_COV'][ref,3,3])
                else :
                  refval = a['PARAM'][ref,6]+refoffset
                  referr = np.sqrt(a['PARAM_COV'][ref,6,6])
                if not usemh: refval -= a['PARAM'][ref,3]
                reflit = (refabun(el,dwarf=dwarf)-refabun('Fe',dwarf=dwarf))
                plots.plotl(ax,xr,[refval-reflit,refval-reflit],color='r')

                # Plot the median of solar abundance stars
                cal=np.where(y[solar] > -9000)[0]
                med = np.median(y[solar[cal]])
                plots.plotl(ax,xr,[med,med],color='y')
                plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.85*(yr[1]-yr[0]),'solar metallicity stars: {:4.2f}'.format(med),color='y')
                if iplot == 0 :
                    out.write(el+'{:8.3f}  {:8d}\n'.format(med,len(cal)))

                # Plot the offset from the optical analysis 
                opt=optabun(el)
                plots.plotl(ax,xr,[opt,opt],color='m')
                plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.75*(yr[1]-yr[0]),'optical offset: {:4.2f}'.format(opt),color='m')

                # Plot M67 points
                clust=fits.open('../../cal/clust.fits')[1].data
                m67=np.where(np.core.defchararray.find(clust['FIELD'],'M67') >= 0)
                m1, m2 = esutil.numpy_util.match(a['APOGEE_ID'][gd],clust['APOGEE_ID'][m67])
                plots.plotp(ax,x[m1],y[m1],size=6,color='k',xr=xr,yr=yr)
                if iplot == 2 :
                    print m1, x[m1]
                if dwarf : 
                    refstar = 'VESTA'
                    if dcal is not None : 
                        refclust=dcal[ielem]  #-dcal[ifeh]
                        plots.plotl(ax,xr,[refclust,refclust],color='g')
                        plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.95*(yr[1]-yr[0]),'M67 dwarfs: {:4.2f}'.format(refclust),color='g')
                else :
                    refstar = 'Arcturus'
                    if gcal is not None : 
                        refclust=gcal[ielem]  #-gcal[ifeh]
                        plots.plotl(ax,xr,[refclust,refclust],color='g')
                        plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.95*(yr[1]-yr[0]),'M67 giants: {:4.2f}'.format(refclust),color='g')
                        if dcal is not None : 
                            drefclust=dcal[ielem] #-dcal[ifeh]
                            plots.plotl(ax,xr,[refclust-drefclust,refclust-drefclust],color='b')
                            plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.90*(yr[1]-yr[0]),'M67 dwarfs: {:4.2f}'.format(drefclust),color='b')
                plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.05*(yr[1]-yr[0]),
                    refstar+'   ASPCAP: {:4.2f}+/-{:4.2f}'.format(refval[0],referr[0])+'   lit: '+'{:4.2f}'.format(reflit),color='r')

            #plt.show()
            plt.savefig(fname+'{:1d}.jpg'.format(iplot))
            plt.close()
            file.append(name+'{:1d}.jpg'.format(iplot))

        ytit.append(name)

#        plt.show()
        files.append(file)
    out.close()
    html.htmltab(files,ytitle=ytit,file='elem/'+prefix+'elem.html',xtitle=xtit,header=tit)

def elemindex() :
    ''' Make the HTML pages (assumes plots have already been made) for individual elements '''

    a,elem,elemtoh = read()

    # loop over elements
    nelem=len(elem)
    for ielem in range(len(elem)+2) :
      if ielem < len(elem) :
          el = elem[ielem].strip()
      elif ielem == nelem :
            el = 'M'
      elif ielem == nelem+1 :
            el = 'alpha'

      ytit=[]
      files=[]
      for prefix in [ 'g','d' ] :
       for suffix in [ '', 'gal' ] :
        name=prefix+el
        file = [prefix+'1mh'+suffix+el+'2.jpg',prefix+'1emh'+suffix+el+'2.jpg',prefix+'1eemh'+suffix+el+'2.jpg',
                prefix+'2mh'+suffix+el+'2.jpg',prefix+'2emh'+suffix+el+'2.jpg',prefix+'2eemh'+suffix+el+'2.jpg',
                prefix+'3mh'+suffix+el+'2.jpg',prefix+'3emh'+suffix+el+'2.jpg',prefix+'3eemh'+suffix+el+'2.jpg']
        xtit = ['Linear 4000-5250','Linear 3750-5250','Linear 3500-5250',
               'Quadratic 4000-5250','Quadratic 3750-5250','Quadratic 3500-5250',
               'Cubic 4000-5250','Cubic 3750-5250','Cubic 3500-5250']
        files.append(file)

      ytit = ['Giants (full)','Giants (70&lt;l&lt;110)','Dwarfs (full)','Dwarfs (70&lt;l&lt;110)']
      html.htmltab(files,file='elem/'+el+'.html',xtitle=xtit,ytitle=ytit)

def main() :
    ''' Make series of plots and web pages for each calibration "type" '''    

    #files = ['testcal','testcal1mh','testcal2mh','testcal3mh',
    #         'testcal1emh','testcal2emh','testcal3emh', 
    #         'testcal1eemh','testcal2eemh','testcal3eemh'] 
    #dirs = ['../testcal','testcal1mh','testcal2mh','testcal3mh',
    #         'testcal1emh','testcal2emh','testcal3emh', 
    #         'testcal1eemh','testcal2eemh','testcal3eemh'] 
    #suffixes = ['','1mh','2mh','3mh','1emh','2emh','3emh','1eemh','2eemh','3eemh']
    files = ['l30e.2']
    dirs = ['../cal']
    suffixes = ['']
    for i in range(len(files)) :
      a,e,etoh = read(file='allStar-'+files[i]+'.fits')
      gcal = fits.open(dirs[i]+'/giantcal.fits')[2].data['ABUN'][0,:,17]
      dcal = fits.open(dirs[i]+'/dwarfcal.fits')[2].data['ABUN'][0,:,17]
      for d in [ False, True ] :
        if d :
          res = ascii.read(dirs[i]+'/dwarfcal.res')
        else :
          res = ascii.read(dirs[i]+'/giantcal.res')
        tmp=etoh   # since etoh gets changed for dwarfs
        plot(a,e,tmp,suffix=suffixes[i],dwarf=d,gcal=gcal,dcal=dcal,res=None,usemh=True,sn=[200,1000])
        plot(a,e,tmp,suffix=suffixes[i]+'gal',dwarf=d,gcal=gcal,dcal=dcal,glon=[70,110],glat=[-5,5],res=None,usemh=True)
    #a,e,etoh = read()
    #plot(a,e,etoh)
    #plot(a,e,etoh,dwarf=True)

def globalscatter(allstar,elems) :
    ''' 
    Compute scatter in clusters
    '''
    clust=apselect.clustdata()
    gd=apselect.select(allstar,badval='STAR_BAD')
    members=[]
    print 'selecting'
    clusts = ['N2420', 'M67', 'N188', 'N7789', 'N6819', 'N6791']
    for cluster in clusts :
        j=np.array(apselect.clustmember(allstar[gd],cluster,raw=True))
        print cluster,len(j)
        members.append(j)
    pdb.set_trace()

    iel=0
    for el in np.append(elems,['M','alpha']) :
        iclust=0
        all=np.array([])
        for cluster in clusts :
            i=np.where(clust.name == cluster)
            mh=clust[i].mh
            name=clust[i].name
            # get cluster members
            j=members[iclust]
            if len(j) > 0 :
                if el.strip() == 'Fe' :
                  abun=allstar['X_H'][gd[j],iel]
                  ok=np.where(((allstar['ELEMFLAG'][gd[j],iel] & 255) == 0) & (allstar['X_H_ERR'][gd[j],iel] < 0.2))[0]
                elif el.strip() == 'M' :
                  abun=allstar['M_H'][gd[j]]
                  ok=np.where(((allstar['PARAMFLAG'][gd[j],3] & 255) == 0) & (allstar['M_H_ERR'][gd[j]] < 0.2))[0]
                elif el.strip() == 'alpha' :
                  abun=allstar['ALPHA_M'][gd[j]]
                  ok=np.where(((allstar['PARAMFLAG'][gd[j],6] & 255) == 0) & (allstar['ALPHA_M_ERR'][gd[j]] < 0.2))[0]
                else :
                  abun=allstar['X_M'][gd[j],iel]
                  ok=np.where(((allstar['ELEMFLAG'][gd[j],iel] & 255) == 0) & (allstar['X_M_ERR'][gd[j],iel] < 0.2))[0]
                if len(ok) > 3 :
                    all=np.append(all,abun[ok]-abun[ok].mean())
            iclust+=1
        print el, all.mean(), all.std(), len(all)
        iel+=1

def getabun(data,elems,el,xh=False,terange=[-1,10000]) :
    '''
    Return the abundance of the requested element, given data array, elem array, element
    '''
    if el.strip() == 'M' :
        ok=np.where(((data['PARAMFLAG'][:,3] & 255) == 0) & (data['FPARAM_COV'][:,3,3] < 0.2) &
                    (data['FPARAM'][:,0] >= terange[0]) & (data['FPARAM'][:,0] <= terange[1]) )[0]
        abun = data['FPARAM'][:,3]
    elif el.strip() == 'alpha' :
        ok=np.where(((data['PARAMFLAG'][:,6] & 255) == 0) & (data['FPARAM_COV'][:,6,6] < 0.2) &
                    (data['FPARAM'][:,0] >= terange[0]) & (data['FPARAM'][:,0] <= terange[1]) )[0]
        abun = data['FPARAM'][:,6]
        if xh : abun+=data['FPARAM'][:,3]
    else :
        iel=np.where(np.core.defchararray.strip(elems['ELEM_SYMBOL'][0]) == el.strip())[0][0]
        ok=np.where(((data['ELEMFLAG'][:,iel] & 255) == 0) & (data['FELEM_ERR'][:,iel] < 0.2) &
                    (data['FPARAM'][:,0] >= terange[0]) & (data['FPARAM'][:,0] <= terange[1]) )[0]
        abun = data['FELEM'][:,iel]
        if xh and not elems['ELEMTOH'][0][iel] : abun+=data['FPARAM'][:,3]
        if not xh and elems['ELEMTOH'][0][iel] : abun-=data['FPARAM'][:,3]
    return abun, ok

def cal(allstar,elems,xh=False,plot=True) :
    ''' 
    Determine internal calibration relations for elements
   
    Args:
        allstar : allStar-like data structure, (i.e., HDU1 of allStar)
        elems : elem-like data structure (e.g. HDU3 of allStar)

    Keyword args:
        xh  : fit in [X/H]? (default=False, i.e. fit in [X/M])
        plot : show individual element plots
    '''

    # select cluster members from array that don't have STAR_BAD into data structure
    clusters=apselect.clustdata()
    clusts = ['M92','M15','M71','N2420', 'M67', 'N188', 'N7789', 'N6819', 'N6791']
    types = [0,1,2,3,4,5,6,7,8]
    markers = ['o','o','o','o','s','s','s','s','s','s']
    colors = ['r','g','b','c','y','m','r','g','b','c']
    gd=apselect.select(allstar,badval='STAR_BAD',raw=True)
    all=[]
    print 'selecting'
    for cluster in clusts :
        j=apselect.clustmember(allstar[gd],cluster,raw=True)
        all=set(all).union(gd[j].tolist())
    data=allstar[list(all)]
    # in the abbreviated array, get the lists of cluster mbmers
    members=[]
    for cluster in clusts :
        j=apselect.clustmember(data,cluster,raw=True)
        members.append(j)
    pdb.set_trace()

    # loop over elements
    rec = np.recarray(len(elems['ELEM_SYMBOL'][0])+2,dtype=[
                       ('elemfit','i4'),
                       ('mhmin','f4'),
                       ('te0','f4'),
                       ('temin','f4'),
                       ('temax','f4'),
                       ('caltemin','f4'),
                       ('caltemax','f4'),
                       ('extfit','i4'),
                       ('extpar','3f4'),
                       ('clust','{:1d}S16'.format(len(clusts))),
                       ('par','3f4'),
                       ('abun','20f4'),
                       ])
    allpars=[]
    iel=0
    if plot :
        fig,ax = plots.multi(1,2,hspace=0.001)
    for el in np.append(elems['ELEM_SYMBOL'][0],['M','alpha']) :
        # parameters for the fit for this element
        pars = dr13cal(el)
        for key in ['elemfit','mhmin','te0','temin','temax','caltemin','caltemax','extfit','extpar'] :
            rec[iel][key]=pars[key]
        rec['clust'] = np.array(clusts,dtype='S16')
        rec['par'] = np.zeros(3)
        if pars['elemfit'] >= 0 :
            # get the good abundance data for this element, load variables for fit (teff, abun, clust)
            abundata, ok = getabun(data,elems,el,xh=xh,terange=[pars['temin'],pars['temax']])
            ind=np.array([],dtype=int)
            clust=np.array([],dtype='S16')
            for iclust in range(len(clusts)) :
                i=np.where(clusters.name == clusts[iclust])
                mh=clusters[i].mh
                if mh > pars['mhmin'] :
                    # get cluster members: intersection of all cluster members and good ones for this element
                    j=list(set(ok).intersection(members[iclust]))
                    if len(j) > 3 :
                        ind=np.append(ind,j)
                        clust=np.append(clust,[clusts[iclust]]*len(j))
            if len(ind) > 0 :
                teff=data['FPARAM'][ind,0]
                abun=abundata[ind]
                visit=data['VISIT'][ind]
                # only use visits=0 for fit
                gd=np.where(visit == 0)[0]
                bd=np.where(visit > 0)[0]
                deriv=calderiv(teff[gd]-pars['te0'],abun[gd],clust[gd],order=pars['elemfit'])
                soln,inv = fit.linear(abun[gd],deriv)
                nclust = len(np.unique(clust[gd]))
                pars['clust'] = np.sort(np.unique(clust[gd]))
                pars['par'] = soln[nclust:len(soln)]
                pars['abun'] = soln[0:nclust]
                func=calfunc(pars,teff,abun,clust,order=pars['elemfit'])
                print '{:<18s}{:8.3f}{:8.3f}'.format(el, (abun[gd]-func[gd]).std(), (abun[bd]-func[bd]).std())
                if plot :
                    ax[0].cla()
                    ax[1].cla()
                    plots.plotp(ax[0],teff[gd],abun[gd]-func[gd], typeref=clust[gd],
                                types=clusts,color=colors,marker=markers,size=16,yt=el)
                    plots.plotp(ax[0],teff[bd],abun[bd]-func[bd],typeref=clust[bd],
                                types=clusts,color=colors,marker=markers,size=16,facecolors='none')
                    func=calfunc(pars,teff,abun,clust,order=0)
                    plots.plotp(ax[1],teff[gd],abun[gd]-func[gd],typeref=clust[gd],
                                types=clusts,color=colors,marker=markers,size=16,xt='Teff',yt=el)
                    plots.plotp(ax[1],teff[bd],abun[bd]-func[bd],typeref=clust[bd],
                                types=clusts,color=colors,marker=markers,size=16,facecolors='none')
                    plots._id_cols=['APOGEE_ID','VISIT']
                    plots._data=data[ind]
                    plots._data_x=teff
                    plots._data_y=abun-func
                    x=np.linspace(pars['temin'],pars['temax'],200)
                    func=calfunc(pars,x,x*0,np.array(['']*len(x)),order=pars['elemfit'])
                    plots.plotl(ax[1],x,func)
                    z=np.where(teff < 4100)[0]
                    for zz in z :
                      print data['APOGEE_ID'][ind[zz]], teff[zz], abun[zz], visit[zz]
                    plots.event(fig)
                    pdb.set_trace()
        allpars.append(pars)

    return allpars

def calfunc(pars,teff,abun,clust,order=1) :
    '''
    Apply calibration function. If clust is not '', then include the mean abundance for the cluster as determined from the fit,
    otherwise only apply the temperature correction

    '''
    npts=len(teff)
    func=np.zeros([npts])
    # if we are given clusters that are not part of the calibration, set them to -999
    j=np.where(clust != '')[0]
    func[j]=-999.
    # start with the cluster mean abundances if requested
    for iclust in range(len(pars['clust'])) :
        j=np.where(clust == pars['clust'][iclust].strip())[0]
        func[j] = pars['abun'][iclust]
    # add the temperature terms, truncating at caltemin and caltemax
    if order >= 1:
        temp=copy.copy(teff)
        bd=np.where(temp < pars['caltemin'])[0]
        temp[bd]=pars['caltemin']
        bd=np.where(temp > pars['caltemax'])[0]
        temp[bd]=pars['caltemax']
        for iorder in range(0,order) :
            func += pars['par'][iorder]*(temp-pars['te0'])**(iorder+1)
    return func

def calderiv(teff,abun,clust,order=1) :
    '''
    Function/derivatives for abundance calibration
    '''
    uclust=np.sort(np.unique(clust))
    npar=order+len(uclust)
    npts=len(teff)
    deriv=np.zeros([npar,npts])
    for iclust in range(len(uclust)) :
        j=np.where(clust == uclust[iclust])[0]
        deriv[iclust,j] = 1.
    if order >= 1:
        for iorder in range(0,order) :
            deriv[len(uclust)+iorder,:] = teff**(iorder+1)
    return deriv
        

def dr13cal(el,dwarfs=False) :
    '''
    Return calibration parameters for requested element

    elemfit gives order/type of polynomial in cluster fit: 1 (linear), 2 (quadratic), 3 (cubic)
    temin/temax gives range over which fit is performed
    caltemin/caltemax gives range over which calibration can be applied (bad outside range)
    extfit gives source of external calibration:  1 (Arcturus), 2 (Vesta), 3 (M67), 4 (solar sequence), 10 (quadratic fit to clusters)
    extpar gives the values of the external calibration
    '''

    # defaults
    te0=4500
    temin=4000
    temax=5000
    elemfit=1
    extfit=0
    caltemin=3532.5
    caltemax=6500
    extpar=[0.,0.,0.]
    mhmin=-1
    if el.strip() == 'Ca' : mhmin = -2.
    if el.strip() == 'C' : mhmin = -0.6
    if el.strip() == 'Fe' : mhmin = -3.
    if el.strip() == 'K' : mhmin = -0.6
    if el.strip() == 'Mn' : mhmin = -2.0
    if el.strip() == 'Na' : mhmin = -0.6
    if el.strip() == 'Ni' : mhmin = -3.0
    if el.strip() == 'N' : mhmin = -0.6
    if el.strip() == 'O' : mhmin = -0.6
    if el.strip() == 'Si' : mhmin = -3.0
    if el.strip() == 'V' : mhmin = -0.6

    if not dwarfs :
        # calibration parameters for giants
        if el.strip() == 'C' :
            elemfit= 0
        elif el.strip() == 'CI' :
            elemfit= 0
        elif el.strip() == 'N' :
            elemfit= 0
        elif el.strip() == 'O' :
            elemfit= 2
            temin= 3750
            extfit= 4
            extpar= [0.060,0.,0.]
        elif el.strip() == 'Na' :
            elemfit= 2
            extfit= 4
            extpar= [0.186,0.,0.]
        elif el.strip() == 'Mg' :
            elemfit= 3
            temin= 3500
            extfit= 4
            extpar= [0.045,0.,0.]
        elif el.strip() == 'Al' :
            elemfit= 3
            extfit= 4
            extpar= [0.108,0.,0.]
        elif el.strip() == 'Si' :
            elemfit= 3
            temin= 3500
            extfit= 4
            extpar= [0.107,0.,0.]
        elif el.strip() == 'P' :
            elemfit= 2
            extfit= 4
            extpar= [-0.008,0.,0.]
        elif el.strip() == 'S' :
            elemfit= 2
            extfit= 4
            extpar= [-0.092,0.,0.]
        elif el.strip() == 'K' :
            elemfit= 1
            extfit= 4
            extpar= [-0.026,0.,0.]
        elif el.strip() == 'Ca' :
            elemfit= 3
            temin= 3750
            extfit= 4
            extpar= [-0.021,0.,0.]
        elif el.strip() == 'Ti' :
            elemfit= 3
            temin= 3500
            extfit= 4
            extpar= [-0.014,0.,0.]
        elif el.strip() == 'TiII' :
            elemfit= 2
            extfit= 4
            extpar= [0.166,0.,0.]
        elif el.strip() == 'V' :
            elemfit= 3
            temin= 3750
            extfit= 4
            extpar= [0.110,0.,0.]
        elif el.strip() == 'Cr' :
            elemfit= 2
            temin= 3500
            extfit= 4
            extpar= [-0.057,0.,0.]
        elif el.strip() == 'Mn' :
            elemfit= 1
            extfit= 4
            extpar= [0.041,0.,0.]
        elif el.strip() == 'Fe' :
            elemfit= 2
            temin= 3500
            extfit= 4
            extpar= [-0.005,0.,0.]
        elif el.strip() == 'Co' :
            elemfit= 3
            extfit= 4
            extpar= [0.003,0.,0.]
        elif el.strip() == 'Ni' :
            elemfit= 2
            temin= 3750
            extfit= 4
            extpar= [-0.001,0.,0.]
        elif el.strip() == 'Cu' :
            elemfit= 3
            temin= 3
            extfit= 4
            extpar= [0.452,0.,0.]
        elif el.strip() == 'Ge' :
            elemfit= 2
            extfit= 4
            extpar= [0.354,0.,0.]
        elif el.strip() == 'Ce' :
            elemfit= -1
        elif el.strip() == 'Rb' :
            elemfit= 2
            temin= 3750
            extfit= 4
            extpar= [-0.105,0.,0.]
        elif el.strip() == 'Y' :
            elemfit= -1
        elif el.strip() == 'Nd' :
            elemfit= -1
        elif el.strip() == 'M' :
            elemfit= 1
        elif el.strip() == 'alpha' :
            elemfit= 2
            extfit= 4
            extpar = [0.056,0.,0.]
    else :

        # default values for dwarfs
        temin=3200
        temax=6250
        elemfit=3
        caltemin=-1
        caltemax=6500
        extfit=0
        extpar=[0.,0.,0.]


        # manual overrides for each element, dwarfs
        if el.strip() == 'C' :
            elemfit=1
            extfit=4
            extpar=[-0.019,0.,0.]
        elif el.strip() == 'CI' :
            extfit=4
            extpar=[-0.026,0.,0.]
        elif el.strip() == 'N' :
            extfit=4
            extpar=[-0.01,0.,0.]
        elif el.strip() == 'O' :
            elemfit=3
            temin=3500
            temax=4500
            extfit=4
            extpar=[0.068,0.,0.]
        elif el.strip() == 'Na' :
            elemfit=1
            temin=3750
            temax=5500
            caltemin=3750
            extfit=4
            extpar=[0.096,0.,0.]
        elif el.strip() == 'Mg' :
            elemfit=3
            temin=3750
            extfit=4
            extpar=[-0.003,0.,0.]
        elif el.strip() == 'Al' :
            elemfit=2
            temin=3750
            caltemin=3500
            extfit=4
            extpar=[0.043,0.,0.]
        elif el.strip() == 'Si' :
            elemfit=1
            temin=3500
            extfit=4
            extpar=[-0.023,0.,0.]
        elif el.strip() == 'P' :
            caltemax=-1
            extfit=0
            extpar=[0.,0.,0.]
        elif el.strip() == 'S' :
            elemfit=1
            temin=3750
            caltemin=5500
            extfit=4
            extpar=[-0.017,0.,0.]
        elif el.strip() == 'K' :
            elemfit=2
            temin=3750
            caltemin=3750
            extfit=4
            extpar=[-0.029,0.,0.]
        elif el.strip() == 'Ca' :
            elemfit=1
            temin=3750
            caltemin=3750
            extfit=4
            extpar=[0.023,0.,0.]
        elif el.strip() == 'Ti' :
            elemfit=3
            temin=3750
            temax=5250
            caltemin=3750
            extfit=4
            extpar=[-0.002,0.,0.]
        elif el.strip() == 'TiII' :
            caltemax=-1
            extfit=0
            extpar=[0.,0.,0.]
        elif el.strip() == 'V' :
            elemfit=2
            temax=5250
            caltemin=3750
            extfit=4
            extpar=[0.002,0.,0.]
        elif el.strip() == 'Cr' :
            elemfit=1
            temax=5250
            caltemin=3750
            extfit=4
            extpar=[-0.044,0.,0.]
        elif el.strip() == 'Mn' :
            elemfit=3
            temin=3500
            caltemin=3500
            extfit=4
            extpar=[-0.077,0.,0.]
        elif el.strip() == 'Fe' :
            elemfit=2
            temin=3500
            extfit=4
            extpar=[0.016,0.,0.]
        elif el.strip() == 'Co' :
            elemfit=-1
        elif el.strip() == 'Ni' :
            elemfit=1
            temin=3500
            caltemin=3500
            extfit=4
            extpar=[0.03,0.,0.]
        elif el.strip() == 'Cu' :
            elemfit=2
            temin=3750
            caltemin=4500
            extfit=4
            extpar=[0.026,0.,0.]
        elif el.strip() == 'Ge' :
            elemfit=-1
        elif el.strip() == 'Ce' :
            elemfit=-1
        elif el.strip() == 'Rb' :
            elemfit=1
            temin=3200
            temax=5250
            extfit=4
            extpar=[-0.217,0.,0.]
        elif el.strip() == 'Y' :
            elemfit=-1
        elif el.strip() == 'Nd' :
            elemfit=-1
        elif el.strip() == 'M' :
            elemfit=3
            temin=3200
            extfit=0
            extpar=[0.0,0.,0.]
        elif el.strip() == 'alpha' :
            elemfit=1
            extfit=4
            extpar=[-0.004,0.,0.]
        
    return {'elemfit': elemfit, 'mhmin' : mhmin, 'te0': te0, 'temin': temin, 'temax': temax, 
            'caltemin': caltemin, 'caltemax' : caltemax, 'extfit' : extfit, 'extpar' : np.array(extpar)}


if __name__ == '__main__' :
    main()
