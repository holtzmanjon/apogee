from holtz.apogee import apselect
from holtz.tools import html
from holtz.tools import match
import os
import shutil
import pdb
import numpy as np
from astropy.io import fits

def calsample(indata,file=None,plot=True,clusters=True,apokasc='APOKASC_cat_v3.6.0',cal1m=True,galcen=True,lowext=True) :
    '''
    selects a calibration subsample from an input apField structure, including several calibration sub-classes: 
        cluster, APOKASC stars, 1m calibration stars. Creates cluster web pages/plots if requested
    '''

    j=np.where(indata['COMMISS'] == 0)[0]
    data=indata[j]
    jc=[]

    if clusters :
        clusts=apselect.clustdata()
        f=html.head(file=file)
        f.write('<TABLE BORDER=2')
        clust=apselect.clustdata()
        for ic in range(len(clust.name)) :
            j=apselect.clustmember(data,clust[ic].name,plot=plot,hard='test')
            print(clust[ic].name,len(j))
            jc.extend(j)
            f.write('<TR><TD>'+clust[ic].name+'<TD>{:12.6f}<TD>{:12.6f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}'.format(
                    clust[ic].ra,clust[ic].dec,clust[ic].rad,clust[ic].rv,clust[ic].drv))
            f.write('<TD><A HREF='+clust[ic].name+'_pos.jpg><IMG SRC='+clust[ic].name+'_pos.jpg width=300></A>')
            f.write('<TD><A HREF='+clust[ic].name+'_pos.jpg><IMG SRC='+clust[ic].name+'_rv.jpg width=300></A>')
            f.write('<TD><A HREF='+clust[ic].name+'_pos.jpg><IMG SRC='+clust[ic].name+'_cmd.jpg width=300></A>')
        html.tail(f)
        print('Number of cluster stars: ',len(jc))
        # remove existing output directories, create new ones
        cleandir('cal/clust',50)
        # create symbolic links in output directories
        nsplit=len(jc)//50+1
        for i in range(len(jc)) :
            symlink(data[jc[i]],'cal/clust',i//nsplit)
        jc=[]
    
    if apokasc is not None :
        apokasc = fits.open(os.environ['IDLWRAP_DIR']+'/data/'+apokasc+'.fits')[1].data
        i1,i2=match.match(data['APOGEE_ID'],apokasc['2MASS_ID'])
        rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
        print('Number of APOKASC RGB stars (every 3rd): ',len(rgb[0:-1:3]))
        jc.extend(i1[rgb][0:-1:3])
        rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
        print('Number of APOKASC RC stars (every 2nd): ',len(rc[0:-2:2]))
        jc.extend(i1[rc][0:-1:2])
        lowg=np.where((apokasc['LOGG_SYD_SCALING'][i2] < 2) & (apokasc['LOGG_SYD_SCALING'][i2] > 0.1))[0]
        print('Number of APOKASC low log g  stars: ',len(lowg))
        jc.extend(i1[lowg])
        highg=np.where((apokasc['LOGG_SYD_SCALING'][i2] > 3.8) & (apokasc['LOGG_SYD_SCALING'][i2] < 5.5))[0]
        print('Number of APOKASC high log g  stars: ',len(highg))
        jc.extend(i1[highg])
        lowz=np.where((apokasc['FE_H_ADOP_COR'][i2] < -1.) & (apokasc['FE_H_ADOP_COR'][i2] > -90.))[0]
        print('Number of APOKASC low [Fe/H] stars: ',len(lowz))
        jc.extend(i1[lowz])
    
    if galcen :
        j=np.where(data['FIELD'] == 'GALCEN')[0]
        print('Number of GALCEN stars: ',len(j))
        jc.extend(j)

    if cal1m :
        j=np.where(data['FIELD'] == 'calibration')[0]
        print('Number of 1m calibration stars: ',len(j))
        jc.extend(j)

    #if lowext :

    # remove existing output directories, create new ones
    cleandir('cal/cal',50)
    # create symbolic links in output directories
    nsplit=len(jc)//50+1
    for i in range(len(jc)) :
        symlink(data[jc[i]],'cal/cal',i//nsplit)

    return jc

def cleandir(out,n) :
    for i in range(n) : 
        try:
            shutil.rmtree('{:s}{:03d}'.format(out,i))
        except : pass
        os.mkdir('{:s}{:03d}'.format(out,i))

def symlink(data,out,idir) :
    outfile='{:s}{:03d}/{:s}.{:04d}.fits'.format(
            out,idir,os.path.splitext(os.path.basename(data['FILE']))[0],data['LOCATION_ID'])
    if data['TELESCOPE'] == 'apo25m' :
        infile='../../{:s}/{:04d}/{:s}'.format(data['TELESCOPE'],data['LOCATION_ID'],data['FILE'])
    else :
        infile='../../{:s}/calibration/{:s}'.format(data['TELESCOPE'],data['FILE'])
    os.symlink(infile,outfile)
