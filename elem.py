# makes a bunch of [X/Fe] vs [Fe/H] plots for use with calibration

import matplotlib.pyplot as plt
import numpy as np
from sdss.apogee import apload
from sdss.apogee import apselect
from holtz.tools import plots
from holtz.tools import html
import pdb
from astropy.io import fits
from astropy.io import ascii
import esutil

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


if __name__ == '__main__' :
    main()
