import os
import numpy as np
import sdss.apogee.apload as apload
import pdb

def location(fields) :
  ''' Return location_ids for requested field(s) '''
  data,hdr=apload.allPlates(hdu=1)
  print('{:<20s}{:<8s}'.format('FIELD','LOCATION_ID'))
  out = ()
  for field in fields :
      print(field)
      j=np.where(np.core.defchararray.strip(data['NAME']) == field.strip())
      print(j)
      n=np.shape(j)[1]
      loc =  int(data['LOCATION_ID'][j[0][0]])
      print('{:<20s}{:<8d}'.format(field,loc))
      out = out + (loc,)
  return out

def plate_mjd(field) :
  ''' List and return plate/mjd combinations for a given field name '''
  data,hdr=apload.allPlates(hdu=1)
  j=np.where(np.core.defchararray.strip(data['NAME']) == field.strip())
  n=np.shape(j)[1]
  print('{:<20s}{:<8s}{:<8s}{:<10s}'.format('FIELD','LOCATION','PLATE','MJD'))
  for i in range(n) :
    print('{:<20s}{:<8d}{:<8d}{:<10d}'.format(data['NAME'][j[0][i]], data['LOCATION_ID'][j[0][i]],data['PLATE'][j[0][i]],data['MJD'][j[0][i]]))
  return data['PLATE'][j],data['MJD'][j]
 
def nums(plate,mjd,type=None) :
  ''' List and return exposure numbers for a given plate/mjd '''
  data,hdr=apload.allExp(hdu=1) 
  if type is None :
      j=np.where((data['PLATEID'] == plate) & (data['MJD'] == mjd) )
  else :
      j=np.where((data['PLATEID'] == plate) & (data['MJD'] == mjd) & (data['IMAGETYP'] == type) )
  print('{:<8s}{:<8s}{:<8s}{:<12s}{:<12s}  {:s}'.format('PLATEID','MJD','CARTID','NUM','DITHPIX','IMAGETYP'))
  n=np.shape(j)[1]
  for i in range(n) :
    print('{:<8d}{:<8d}{:<8d}{:<12d}{:<12.3f}  {:s}'.format(data['PLATEID'][j[0][i]], data['MJD'][j[0][i]],data['CARTID'][j[0][i]],data['NUM'][j[0][i]],data['DITHPIX'][j[0][i]],data['IMAGETYP'][j[0][i]]))
  return data['NUM'][j]
