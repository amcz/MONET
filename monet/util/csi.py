#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *
import sys 
from scipy.io import netcdf
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import string
import datetime
#import shapely.geometry as sgeo
from hysplit import *
import itertools

"""
routines to help calculates critical success index by a point by point comparison.
data must be gridded the same.

Functions
---------
match_ra
get_area_ra
calc_fss  (fraction skill score)
find_threshold
mask_threshold
calc_csi (critical success index, aka figure of merit in space)



"""
def match_ra(ra1, lat1, lon1, ra2, lat2, lon2, missing_value=0, verbose=1):
    """inputs two data 'structures' data, latitude, longitude.
    ra1 : numpy array 
    ra2 : numpy array 
    lon1 : numpy array 
    lat1 : numpy array 
    lon2 : numpy array 
    lat2 : numpy array 
    shape of ra1, lat1, lon1 should be the same
    shape of ra2, lat2, lon2 should be the same.
    Assumes lat-lon grids are spaced the same but may be offset. 
    output arrays can be input into calc_csi 

    return
    latlongrid 
    ra1 
    ra2
    """
    latlongrid=[]
    if verbose==1: 
        print 'RA1 ' , ra1.shape , '\n LAT1 ' , lat1.shape , '\n LON1 ', lon1.shape
        print 'RA2 ' , ra2.shape , '\n LAT2 ' , lat2.shape , '\n LON2 ', lon2.shape
    ##check to make shapes of input arrays are the same.
    if ra1.shape == lat1.shape and ra1.shape == lat1.shape and  \
       ra2.shape == lat2.shape and ra2.shape == lat2.shape:
        
       ##check to make sure data array not already equal. 
       if not (np.array_equal(lon1, lon2) and np.array_equal(lat1, lat2)):
          
           ###MATCH LATITUDES######## 
           lat_space =  round(lat1[1][0] - lat1[0][0], 3)
           lat_space2 = round(lat2[1][0] - lat2[0][0], 3)
           if verbose ==1:
              print 'Lat spacing', lat_space , lat_space2
           #print lat1 
           if lat_space == 0:
              print 'WARNING. latitude spacing is zero!'
           if lat_space != lat_space2:
              print 'WARNING: spacing between two arrays is not the same'
           minlat = np.amin(np.append(lat1, lat2))
           maxlat = np.amax(np.append(lat1, lat2))
           if verbose == 1:
               print 'Min latitude ' , minlat , ' ' , np.amin(lat1),  ' ' , np.amin(lat2) 
               print 'Max latitude ' ,  maxlat, ' ' , np.amax(lat1),  ' ' , np.amax(lat2)
           latlist = np.arange(minlat, maxlat+lat_space, lat_space)
 

           ##Add padding to ra1 ########
           la1 = int(round((np.amin(lat1)-np.amin(latlist))/lat_space))
           la2 = int(round((np.amax(latlist)-np.amax(lat1))/lat_space))
           tempra = np.zeros(len(lat1[0,:])).reshape(1,len(lat1[0,:]))
           tempra = tempra + missing_value
           for x in range(la2):                                    
               ra1 = np.concatenate((ra1, tempra), axis=0)        
           for x in range(la1):
               ra1 = np.concatenate((tempra, ra1), axis=0)
           #print 'RA1 \n' , ra1 
 
           ##Add padding to ra2 ########
           la1 = int(round((np.amin(lat2)-np.amin(latlist))/lat_space))
           la2 = int(round((np.amax(latlist)-np.amax(lat2))/lat_space))
           #print la1, la2
           tempra = np.zeros(len(lat2[0,:])).reshape(1,len(lat2[0,:]))
           tempra = tempra + missing_value
           for x in range(la2): 
               if verbose:
                  print 'la2', x                                   
               ra2 = np.concatenate((ra2, tempra), axis=0)        
           for x in range(la1):
               if verbose:
                  print 'la1' , x                                   
               ra2 = np.concatenate((tempra, ra2), axis=0)
           #print 'RA2 \n' , ra2 


           ###MATCH LONGITUDES######## 
           lon_space =  round(lon1[0][1] - lon1[0][0] , 3)
           lon_space2 = round(lon2[0][1] - lon2[0][0] , 3)
           if verbose ==1:
              print 'LON space' , lon_space , lon_space2
           if lon_space == 0:
              print 'WARNING. longitude spacing is zero!'
           if lon_space != lon_space2:
              print 'WARNING: longitude spacing for two arrays is not the same.'
           minlon = np.amin(np.append(lon1, lon2))
           maxlon = np.amax(np.append(lon1, lon2))
           lonlist = np.arange(minlon, maxlon+lon_space, lon_space)
           if verbose ==1:
              print 'Min longitude ' , minlon , 'Max longitude ', maxlon
      
           ##add padding to ra1 ################
           la1 = int(round((np.amin(lon1)-np.amin(lonlist))/lon_space))
           la2 = int(round((np.amax(lonlist)-np.amax(lon1))/lon_space))
           #print la1, la2
           tempra = np.zeros(len(ra1[:,0])).reshape(1,len(ra1[:,0]))
           tempra = tempra + missing_value
           #print 'ra shape'   , ra1.shape
           for x in range(la2):
               ra1=np.concatenate((ra1,tempra.T), axis=1)
           for x in range(la1):
               ra1=np.concatenate((tempra.T, ra1), axis=1)
           #print 'RA1 \n' , ra1
        
           ##add padding to ra2 ################
           la1 = int(round((np.amin(lon2)-np.amin(lonlist))/lon_space))
           la2 = int(round((np.amax(lonlist)-np.amax(lon2))/lon_space))
           #print la1, la2
           tempra = np.zeros(len(ra2[:,0])).reshape(1,len(ra2[:,0]))
           tempra = tempra + missing_value
           #print 'ra shape'   , ra2.shape
           for x in range(la2):
               ra2=np.concatenate((ra2,tempra.T), axis=1)
           for x in range(la1):
               ra2=np.concatenate((tempra.T, ra2), axis=1)
           #print 'RA2 \n' , ra2
           latlongrid = np.meshgrid(lonlist, latlist)
       else:
           latlongrid = np.meshgrid(lon1, lat1)    
           #print latlongrid[0]
           #print latlongrid[1]
 
    else:
       print "Error: shapes of input arrays are incorrect."
    if verbose ==1:
       print 'Shapes of output arrays ra1, ra2, lat, lon' , ra1.shape, ra2.shape, latlongrid[0].shape, latlongrid[1].shape
    return latlongrid , ra1 , ra2


def get_area_ra(lat, lon, radius=6378.137):
    """Input is a 2d latitude numpy array and a 2d longitude numpy array. 
       Output is array of same size with corresponding area (km2) of each grid cell.
       Converts degrees to meters using a radius of 6378.137 km. Can  use the alt keyword to change this radius.
       Assumes grid is spaced equally in degrees lat and degrees lon"""

    d2km = radius * pi / 180.0     #convert degree latitude to meters. Assumes Earth radius of 6378137 meters.
    d2r =  pi/180.0             #convert degrees to radians

    lat_space =  lat[1][0] - lat[0][0]
    lon_space =  lon[0][1] - lon[0][0]
    if lat_space == 0 or lon_space==0:
       print 'WARNING: lat or lon spacing is zero'
    area = lat_space*d2km * lon_space * d2km * cos(d2r * lat)
    #print 'Area shape', area.shape , lat.shape , np.amin(area) , np.amax(area) 
    return area


def calc_fss(ra1, ra2, nodata_value='', threshold=0, verbose=0, sn=[0,3]):
   """calculates the fss (fraction skill score See Robers and Lean, Monthly
      Weather Review, V126 2007..)

       sn is the size of the squares to use in number of pixels to a side.
       default is to use 1, 3, 5, 7, 9 pixels size squares
   """
   fss_list = []
   print 'RANGE' , range(sn[0],sn[1])
   bigN = ra1.size
   for sz in range(sn[0],sn[1]):
       ijrange = range(-1*sz, sz+1)
       x = itertools.permutations(ijrange+ijrange,2)
       li=[]
       for i  in x:
           li.append(i)
       incrlist=list(set(li))
       square_size = len(incrlist)
       print 'square size ' , sz*2+1 , square_size
       icol = ra1.shape[1]
       irow = ra1.shape[0]
       maxrow = irow
       maxcol = icol
       print ra1.shape , irow , icol
       mra1 = mask_threshold(ra1, threshold=0)
       mra2 = mask_threshold(ra2, threshold=threshold)
       print 'SUMS' , mra1.sum() , mra2.sum()
       list_fractions_1 = np.empty_like(mra1).astype(float)
       list_fractions_2 = np.empty_like(mra1).astype(float)
       for ic in range(0,icol):
           for ir in range(0,irow):
               fraction_1=0
               fraction_2=0
               for incr in incrlist:
                   if ir+incr[0] < maxrow and ic+incr[1] < maxcol and ir+incr[0]>=0 and ic+incr[1] >=0:
                      fraction_1 += mra1[ir + incr[0]][ic+incr[1]]
                      fraction_2 += mra2[ir + incr[0]][ic+incr[1]]
               fraction_1 = float(fraction_1) / float(square_size)          
               fraction_2 = float(fraction_2) / float(square_size)          
               list_fractions_1[ir][ic] = float(fraction_1)
               list_fractions_2[ir][ic] = float(fraction_2)
       fbs = np.power(list_fractions_1 - list_fractions_2, 2).sum() / float(bigN)
       print 'FBS ' , fbs
       fbs_ref = (np.power(list_fractions_1,2).sum() + np.power(list_fractions_2,2).sum() ) / float(bigN)
       print 'FBS reference' , fbs_ref
       fss = 1 - (fbs / fbs_ref)
       print 'FSS ' , fss
       fss_list.append((fss, len(incrlist)))
   return fss_list


def find_threshold(ra1, ra2, nodata_value=None):
    """ra1 is the satellite data.
       ra2 is model data"""
    mask1 = mask_threshold(ra1, threshold=0)
    #numpixels = mask1.size - mask1.sum()
    numpixels = mask1.sum()
    list2 = np.copy(ra2)
    if nodata_value:
       vpND = np.where(ra1 == nodata_value)
       #print 'vpND' , len(vpND[0]) , vpND
       vpNDB = np.where(list2[vpND] > 0)
       #print len(vpNDB[0])
       list2[vpND] = 0 

    list2 = list2.reshape(list2.shape[0] * list2.shape[1]) 
    list2.sort()
    list2 = list2[::-1]
    np.set_printoptions()
    print 'match treshold' , numpixels, list2[numpixels] , list2[numpixels + 1000] , list2[numpixels-1000]
    print mask1.size , mask1.sum() , numpixels 
    mask2 = mask_threshold(list2, threshold=0)
    vpi = np.where(list2 > 0)
    print list2.size  , mask2.sum() , np.min(list2[vpi])

    return list2[numpixels]

def mask_threshold(ra1, threshold=0):
    """input array and threshold. Returns array with 1's where value of array above threshold and 0 otherwise"""
    vp1 = np.where(ra1 > threshold)
    mask1 = np.zeros_like(ra1)
    mask1[vp1] = 1
    return mask1

def calc_csi(ra1, ra2, area=[], nodata_value='', threshold=0, verbose=0):
    """ra1 and ra2 and area need to have the same shape. Assumes each point in array at same position.
       If the points have different areas associated with them then an area array can be input as well.
       ra1 should be the observations/satellite. ra2 should be the forecast/model."""
    area=np.array(area)
    csihash = {}
    #vp1 = np.where(ra1 > threshold)
    vp1 = np.where(ra1 > 0)
    vptest = np.where(logical_and(ra1>0 , ra1 <= threshold))
    print 'TEST VPTEST ' , threshold , len(vptest)
    #print 'VP1' , vp1
    mask1 = ra1[:] * 0
    mask1[vp1] = 1
    if threshold != 0:
       vp2 = np.where(ra2 >= threshold)
    else:
       vp2 = np.where(ra2 > threshold)
    mask2 = ra2[:] * 0
    mask2[vp2] = 1

    matchra = mask2 * mask1  #returns array with ones where the two arrays both have valid values.
    onlyra1 = mask1 - matchra  #returns array with ones where ra1 has valid values and ra2 does not.
    onlyra2 = mask2 - matchra  #returns array with ones where ra2 has valid values and ra1 does not.
    csihash['matched'] = matchra
    csihash['ra1'] = onlyra1
    csihash['ra2'] = onlyra2
    allra = matchra + onlyra1 + onlyra2 ##ra with 1's in union. 0's where both are 0.
    vpi = np.where(allra==1)
    ##Find pattern correlation (from Zidikheri and Potts) doesn't make sense to me.
    totalpts = matchra.sum() + onlyra1.sum() + onlyra2.sum()
    totalpts = mask2.shape[0] * mask2.shape[1] 
    ra1ave = (onlyra1.sum() + matchra.sum() ) / float(totalpts)
    ra2ave = (onlyra2.sum() + matchra.sum() ) / float(totalpts)
    #ra1corr = (mask1 - ra1ave) * allra
    #ra2corr = (mask2 - ra2ave) * allra
    ra1corr = (mask1 - ra1ave) 
    ra2corr = (mask2 - ra2ave)
    norm =((ra1corr *ra1corr).sum())**0.5 * ((ra2corr*ra2corr).sum())**0.5 
    pcorr = (ra1corr * ra2corr).sum()  / norm
    print 'PCORR' , pcorr 
    print 'ra1ave (obs)' ,  ra1ave 
    print 'ra2ave (calc)' , ra2ave, 
    print 'NORM' , norm , 'ra1corr*ra2corr' ,  (ra1corr *ra2corr).sum()
#    print totalpts, mask1.shape
    print mask1.sum(), mask2.sum(), np.max(ra1corr), np.min(ra1corr)

    ra1ave = 0
    ra2ave = 0
    ra1corr = (mask1 - ra1ave) 
    ra2corr = (mask2 - ra2ave)
    norm =((ra1corr *ra1corr).sum())**0.5 * ((ra2corr*ra2corr).sum())**0.5 
    pcorr = (ra1corr * ra2corr).sum()  / norm
    print 'PCORR (uncentered)' , pcorr
    #plt.imshow(mask1)
    #plt.show()
    #plt.imshow(mask2)
    #plt.show()

    #pcorr=-999 
    ##Find where data arrays have no information.
    if nodata_value != '':
       vpND = np.where(logical_or(ra1 == nodata_value , ra2==nodata_value))
       maskND = ra2[:] * 0
       maskND[vpND] = 1
       vp_ra1ND = np.where(logical_and(maskND==1 , onlyra1==1))     
       if vp_ra1ND[0] != []:  
          onlyra1[vp_ra1ND] = 0
       vp_ra2ND = np.where(logical_and(maskND==1 , onlyra2==1))     
       if vp_ra2ND[0] != []:  
          onlyra2[vp_ra2ND] = 0

    if area.shape == matchra.shape:
       csihash['CSI'] = (matchra*area).sum() / ((matchra*area).sum() + (onlyra1*area).sum() + (onlyra2*area).sum())
       csihash['POD'] = (matchra*area).sum() / ((matchra*area).sum() + (onlyra1*area).sum())
       csihash['FAR'] = (onlyra2*area).sum() / ((matchra*area).sum() + (onlyra2*area).sum())
       if verbose ==1:
         print 'used area'
         print (matchra*area).sum() , (onlyra1*area).sum() , (onlyra2*area).sum()
         print  'CSI POD FAR' , csihash['CSI'] , csihash['POD'] , csihash['FAR']
    else: 
       csihash['CSI'] = matchra.sum() / (matchra.sum() + onlyra1.sum() + onlyra2.sum())
       csihash['POD'] = matchra.sum() / (matchra.sum() + onlyra1.sum())
       csihash['FAR'] = onlyra2.sum() / (matchra.sum() + onlyra2.sum())
       csihash['PCORR'] = pcorr
       if verbose ==1:
          print 'HERE' , matchra.size, matchra.shape , onlyra2.shape
          print  matchra.sum() , onlyra1.sum() , onlyra2.sum()
          print  'CSI POD FAR' , csihash['CSI'] , csihash['POD'] , csihash['FAR']
 
    return csihash

    






