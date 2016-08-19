import matplotlib.pyplot as plt
import numpy as np
import pdb
import sdss.apogee.apload as apload

def plot1d(hd,row,erase=True,col=False,chip=True,pixel=False) :
  """
    NAME: applot.plot1d
    PURPOSE:  Plot a single row of input HDUList
    USAGE:  applot.plot1d(HDUList,row[,erase=False,col=True,chip=False])
  """
  if erase :
    plt.clf()

  colors=['r','g','b']
  i=0
  for chip in ['a','b','c'] :
    if col :
      plt.plot(hd[chip][1].data[:,row],colors[i])
    else :
      if pixel :
          plt.plot(np.arange(2048)+i*2048,hd[chip][1].data[row,:],colors[i])
      else :
          plt.plot(hd[chip][4].data[row,:],hd[chip][1].data[row,:],colors[i])
    i+=1

  plt.draw()

def apstar(hd,row=0,erase=True) :
  if erase :
    plt.clf()
  wave=hd[1].header['CRVAL1']+np.arange(hd[1].header['NAXIS1'])*hd[1].header['CDELT1']
  wave=10**wave
  plt.plot(wave,hd[1].data[row,:])

def aspcapStar(hd,erase=True) :
  if erase :
    plt.clf()
  wave=hd[1].header['CRVAL1']+np.arange(hd[1].header['NAXIS1'])*hd[1].header['CDELT1']
  wave=10**wave
  plt.plot(wave,hd[1].data)
  plt.plot(wave,hd[3].data)
  plt.plot(wave,hd[2].data,'r')
  plt.ylim(0,1.3)
  plt.draw()

def compare(hd1,hd2,row,v5=None,v6=None,chip=True) :
  """
    NAME: applot.compare
    PURPOSE:  Do a comparison plot of 2 different HDUList 
    USAGE:  applot.compare(HDUList1,HDUList2,row[,chip=False])
  """
  ax1 = plt.subplot2grid((3,1), (0,0))
  ax1.set_ylim([0.9,1.1])
  ax1.set_xlabel('Wavelength')
  ax1.set_ylabel('Flux ratio')
  ax2 = plt.subplot2grid((3,1), (1,0), rowspan=2)
  ax2.set_xlabel('Wavelength')
  ax2.set_ylabel('Flux')
  if chip :
    med=np.median(hd1['b'][1].data[row,:])
  else :
    med=np.median(hd1[1].data)
  if med < 100 :
    med = 100.
  ax2.set_ylim(-50.,2.0*med)

  colors=['r','g','b']
  i=0
  if chip :
    for chip in ['a','b','c'] :
      ax1.plot(hd1[chip][4].data[row,:],hd1[chip][1].data[row,:]/hd2[chip][1].data[row,:],colors[i]) 
      ax1.plot(hd1[chip][4].data[row,:],hd1[chip][5].data[row,:],'g') 
      ax1.plot(hd1[chip][4].data[row,:],hd1[chip][7].data[row,:],'b',linestyle='dotted')
      ax2.plot(hd1[chip][4].data[row,:],hd1[chip][1].data[row,:],'r')
      ax2.plot(hd2[chip][4].data[row,:],hd2[chip][1].data[row,:],'g')
      if v5 is None :
          ax2.text(15100,0.2*med,hd1[chip][11].data['OBJECT'][row])
      else :
          j=np.where(v5['FIBERID'] == hd1[chip][11].data['FIBERID'][row])
          if len(j[0]) > 0 :
            ax2.text(15100,0.2*med,str(hd1[chip][11].data['OBJECT'][row])+'  RVs: {:7.2f}{:7.2f}'.format(v5['VREL'][j[0][0]],v6['VREL'][j[0][0]]))
          else :
            ax2.text(15100,0.2*med,hd1[chip][11].data['OBJECT'][row])
  else :
      for i in range(3) :
        ax1.plot(hd1[4].data[i,:],hd1[1].data[i,:]/hd2[1].data[i,:],colors[i]) 
        ax1.plot(hd1[4].data[i,:],hd1[5].data[i,:],'g') 
        ax1.plot(hd1[4].data[i,:],hd1[7].data[i,:],'b',linestyle='dotted')
        ax2.plot(hd1[4].data[i,:],hd1[1].data[i,:],'r')
        ax2.plot(hd1[4].data[i,:],hd2[1].data[i,:],'g')

def platecomp(plate,mjd,vers1,vers2) :
    """
      NAME: applot.platecomp
      PURPOSE:  Do a comparison plot of 2 different HDUList 
      USAGE:  applot.platecomp(plate,mjd,vers1,vers2)
    """
    apload.apred = vers1
    r5=apload.apPlate(plate,mjd)
    location=r5['a'][0].header['LOCID']
    v5,h=apload.apVisitSum(location,plate,mjd,hdu=1)
    apload.apred = vers2
    r6=apload.apPlate(plate,mjd)
    location=r6['a'][0].header['LOCID']
    v6,h=apload.apVisitSum(location,plate,mjd,hdu=1)
    for i in range(300) :
      compare(r5,r6,299-i,v5,v6)
      print('Use matplotlib window tools to zoom if desired')
      print('Hit c to continue to next spectrum')
      plt.draw()
      pdb.set_trace()

def comp1m(program,mjd,obj,vers1,vers2) :
    ''' Compare 1m spectrum for a given PROGRAM/MJD/OBJECT for two input versions '''
    apload.apred = vers1
    r5=apload.apVisit1m(program,mjd,obj)
    apload.apred = vers2
    r6=apload.apVisit1m(program,mjd,obj)
    for i in range(300) :
      compare(r5,r6,299-i,chip=False)
      pdb.set_trace()
