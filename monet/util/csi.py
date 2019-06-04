#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *
import sys 
#from scipy.io import netcdf
#from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import string
import datetime
#import shapely.geometry as sgeo
#from hysplit import *
from itertools import permutations
import time


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
        print('RA1 ' , ra1.shape , '\n LAT1 ' , lat1.shape , '\n LON1 ', lon1.shape)
        print('RA2 ' , ra2.shape , '\n LAT2 ' , lat2.shape , '\n LON2 ', lon2.shape)
    ##check to make shapes of input arrays are the same.
    if ra1.shape == lat1.shape and ra1.shape == lat1.shape and  \
       ra2.shape == lat2.shape and ra2.shape == lat2.shape:
        
       ##check to make sure data array not already equal. 
       if not (np.array_equal(lon1, lon2) and np.array_equal(lat1, lat2)):
          
           ###MATCH LATITUDES######## 
           lat_space =  round(lat1[1][0] - lat1[0][0], 3)
           lat_space2 = round(lat2[1][0] - lat2[0][0], 3)
           if verbose ==1:
              print('Lat spacing', lat_space , lat_space2)
           #print lat1 
           if lat_space == 0:
              print('WARNING. latitude spacing is zero!')
           if lat_space != lat_space2:
              print('WARNING: spacing between two arrays is not the same')
           minlat = np.amin(np.append(lat1, lat2))
           maxlat = np.amax(np.append(lat1, lat2))
           if verbose == 1:
               print('Min latitude ' , minlat , ' ' , np.amin(lat1),  ' ' , np.amin(lat2)) 
               print('Max latitude ' ,  maxlat, ' ' , np.amax(lat1),  ' ' , np.amax(lat2))
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
                  print('la2', x)                                   
               ra2 = np.concatenate((ra2, tempra), axis=0)        
           for x in range(la1):
               if verbose:
                  print('la1' , x)                                   
               ra2 = np.concatenate((tempra, ra2), axis=0)
           #print 'RA2 \n' , ra2 


           ###MATCH LONGITUDES######## 
           lon_space =  round(lon1[0][1] - lon1[0][0] , 3)
           lon_space2 = round(lon2[0][1] - lon2[0][0] , 3)
           if verbose ==1:
              print('LON space' , lon_space , lon_space2)
           if lon_space == 0:
              print('WARNING. longitude spacing is zero!')
           if lon_space != lon_space2:
              print('WARNING: longitude spacing for two arrays is not the same.')
           minlon = np.amin(np.append(lon1, lon2))
           maxlon = np.amax(np.append(lon1, lon2))
           lonlist = np.arange(minlon, maxlon+lon_space, lon_space)
           if verbose ==1:
              print('Min longitude ' , minlon , 'Max longitude ', maxlon)
      
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
       print("Error: shapes of input arrays are incorrect.")
    if verbose ==1:
       print('Shapes of output arrays ra1, ra2, lat, lon' , ra1.shape, ra2.shape, latlongrid[0].shape, latlongrid[1].shape)
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
       print('WARNING: lat or lon spacing is zero')
    area = lat_space*d2km * lon_space * d2km * cos(d2r * lat)
    #print 'Area shape', area.shape , lat.shape , np.amin(area) , np.amax(area) 
    return area


def calc_fss(ra1, ra2, nodata_value='', threshold=0, verbose=0, sn=[0,5]):
   """Calculates the fraction skill score (fss) 
       See Robers and Lean (2008) Monthly Weather Review
       and Schwartz et al (2010) Weather and Forecasting
       for more information.
       
       ra1 is observations/satellite
       ra2 is the forecast/model
       threshold = value for data threshold
       sn is the number of pixels (radius) to use in fractions calculation
            default is to use 1, 3, 5, 7, 9 pixels size squares
   """
   fss_list = []
   if (verbose == 1):
       print('RANGE' , list(range(sn[0],sn[1])))
   bigN = ra1.size
   # create binary fields
   mra1 = mask_threshold(ra1, threshold = 0)
   mra2 = mask_threshold(ra2, threshold = threshold)
   for sz in range(sn[0],sn[1]):
       # e.g. for sz=3, ijrange=[-3,-2,-1,0,1,2,3]
       ijrange = list(range(-1*sz, sz+1))
       print('ijrange: ', ijrange)
       # li will be coordinates of all the points in the ra.
       # e.g. [(-3,-3), (-3,-2), (-3,-1).....]
       x = permutations(ijrange+ijrange,2)
       li = []
       for i in x:
           li.append(i)
       # remove duplicates (e.g. (2,3) will be in there twice)
       incrlist=list(set(li))
       print(incrlist)
       square_size = len(incrlist)
       if (verbose == 1):
           print('square size ' , sz*2+1 , square_size)
           print(ra1.shape)
       
       # total number of above threshold points.
       print('SUMS' , mra1.sum() , mra2.sum())
       icol = ra1.shape[1]
       irow = ra1.shape[0]
       maxrow = irow
       maxcol = icol
       # create empty arrays of same shape as mra1 and mra2
       list_fractions_1 = np.empty_like(mra1).astype(float)
       list_fractions_2 = np.empty_like(mra1).astype(float)
       start = time.time()
       for ic in range(0,icol):
           for ir in range(0,irow):
               fraction_1 = 0
               fraction_2 = 0
               for incr in incrlist:
                   #Finding pixels for fraction calculation
                   #Requires designated box for calculation to be fully within domain
                   #doesn't calculate at edges
                   if ir+incr[0] < maxrow and ic+incr[1] < maxcol and ir+incr[0] >= 0 and ic+incr[1] >= 0:
                       fraction_1 += mra1[ir + incr[0]][ic+incr[1]]
                       fraction_2 += mra2[ir + incr[0]][ic+incr[1]]
               #if fraction_1 > 0 or fraction_2 > 0:
               list_fractions_1[ir][ic] = (float(fraction_1) / float(square_size))
               list_fractions_2[ir][ic] = (float(fraction_2) / float(square_size))
       end = time.time()
       print('Total time: ', end - start)
       #Can plot fractions if desired (double check calculations)
       plot_fractions = False
       if (plot_fractions == True):
          fig = plt.figure(1)
          ax1 = fig.add_subplot(2,1,1) 
          ax2 = fig.add_subplot(2,1,2) 
          ax1.imshow(list_fractions_1)
          ax2.imshow(list_fractions_2)
          plt.show()
       #Calculate the Fractions Brier Score (FBS)
       fbs = np.power(list_fractions_1 - list_fractions_2, 2).sum() / float(bigN)
       print('FBS ' , fbs)
       #Calculate the worst possible FBS (assuming no overlap of nonzero fractions)
       fbs_ref = (np.power(list_fractions_1,2).sum() + np.power(list_fractions_2,2).sum() ) / float(bigN)
       print('FBS reference' , fbs_ref)
       #Calculate the Fractional Skill Score (FSS)
       fss = 1 - (fbs / fbs_ref)
       print('FSS ' , fss)
       fss_list.append((fss, len(incrlist)))
       print(' ')
   return fss_list

def find_threshold(ra1, ra2, nodata_value=None):
    """
       Base threshold on matching number of pixels in observations.
       ra1 is the satellite data.
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
    print('match treshold' , numpixels, list2[numpixels] , list2[numpixels + 1000] , list2[numpixels-1000])
    print(mask1.size , mask1.sum() , numpixels) 
    mask2 = mask_threshold(list2, threshold=0)
    vpi = np.where(list2 > 0)
    print(list2.size  , mask2.sum() , np.min(list2[vpi]))
    return list2[numpixels]

def mask_threshold(ra1, threshold = 0):
    """input array and threshold. 
        Returns array with 1's where value is above threshold, and 0 otherwise"""
    vp1 = np.where(ra1 > threshold)
    mask1 = np.zeros_like(ra1)
    mask1[vp1] = 1
    return mask1

def trim_arrays(mask1, mask2, sn=[0,5]):
    """Trims binary mask arrays to array size that encompasses data and minor
    padding for both mask input arrays. Requires pixel radius (sn) value to determine 
    padding to keep in array. Will reduce time of fss calculation."""

    pad = sn[1]
    rows1 = np.any(mask1, axis=1)
    cols1 = np.any(mask1, axis=0)
    rows2 = np.any(mask2, axis=1)
    cols2 = np.any(mask2, axis=0)
    ymin1, ymax1 = np.where(rows1)[0][[0, -1]]
    xmin1, xmax1 = np.where(cols1)[0][[0, -1]]
    ymin2, ymax2 = np.where(rows2)[0][[0, -1]]
    xmin2, xmax2 = np.where(cols2)[0][[0, -1]]
    trim1 = mask1[(ymin1-pad):(ymax1+1+pad), (xmin1-pad):(xmax1+1+pad)]
    trim2 = mask2[(ymin2-pad):(ymax2+1+pad), (xmin2-pad):(xmax2+1+pad)]

    yarr = [ymin1, ymin2, ymax1, ymax2]
    xarr = [xmin1, xmin2, xmax1, xmax2]
    trimboth1 = mask1[(min(yarr)-pad):(max(yarr)+1+pad), (min(xarr)-pad):(max(xarr)+1+pad)]
    trimboth2 = mask2[(min(yarr)-pad):(max(yarr)+1+pad), (min(xarr)-pad):(max(xarr)+1+pad)]

    return trim1, trim2, trimboth1, trimboth2

def calc_csi(ra1, ra2, area=None, nodata_value='', threshold=0, verbose=0):
    """ra1 and ra2 need to be on the same grid. 
        See monet.remap_nearest or monet.remap_xesmf to remap arrays.
        ra1 = observations/satellite
        ra2 = the forecast/model
        area = optional array of grid areas
        nodata_value = flag for expanded array grid cells created if ash near boundary
        threshold = value for data threshold
        CSI equation: hits / (hits + misses + false alarms)"""
    #ra1 = ra1.load()
    #ra2 = ra2.load()
    area=np.array(area)
    #Creating a csihash disctionary
    csihash = {}
    #vptest = np.where(logical_and(ra1 > 0 , ra1 <= threshold))
    #if (verbose == 1):
    #    print('TEST VPTEST ' , threshold , vptest)

    #Converting ra1 and ra2 to arrays of 0's and 1's (1 with values, 0 no values)
    mask1 = mask_threshold(ra1, threshold = 0)
    if (threshold != 0):
       vp2 = np.where(ra2 >= threshold)
       mask2 = ra2[:] * 0
       mask2[vp2] = 1
    else:
       mask2 = mask_threshold(ra2, threshold = 0)

    #Calculating hits (matchra), misses (onlyra1), false alarms (onlyra2) for CSI calculation
    matchra = mask2 * mask1  #returns array with ones where the two arrays both have valid values.
    onlyra1 = mask1 - matchra  #returns array with ones where ra1 has valid values and ra2 does not.
    onlyra2 = mask2 - matchra  #returns array with ones where ra2 has valid values and ra1 does not.
    #Assigning a, b, and c arrays to csihash dictionary
    csihash['matched'] = matchra
    csihash['ra1'] = onlyra1
    csihash['ra2'] = onlyra2

    allra = matchra + onlyra1 + onlyra2 ##ra with 1's in union. 0's where both are 0. (THIS ISNT TRUE)
    #allra would have 3 at union, 1 at onlyra1, 1 at onlyra2, and 0 where both are 0, right???
    vpi = np.where(allra == 1) ## Why is this necessary?
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
    print('PCORR' , pcorr) 
    print('ra1ave (obs)' ,  ra1ave) 
    print('ra2ave (calc)' , ra2ave, end=' ') 
    print('NORM' , norm , 'ra1corr*ra2corr' ,  (ra1corr *ra2corr).sum())
#    print totalpts, mask1.shape
    print(mask1.sum(), mask2.sum(), np.max(ra1corr), np.min(ra1corr))

    ra1ave = 0
    ra2ave = 0
    ra1corr = (mask1 - ra1ave) 
    ra2corr = (mask2 - ra2ave)
    norm =((ra1corr *ra1corr).sum())**0.5 * ((ra2corr*ra2corr).sum())**0.5 
    pcorr = (ra1corr * ra2corr).sum()  / norm
    print('PCORR (uncentered)' , pcorr)
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
         print('used area')
         print((matchra*area).sum() , (onlyra1*area).sum() , (onlyra2*area).sum())
         print('CSI POD FAR' , csihash['CSI'] , csihash['POD'] , csihash['FAR'])
    else: 
       csihash['CSI'] = matchra.sum() / (matchra.sum() + onlyra1.sum() + onlyra2.sum())
       # hit rate or probability of detection (p 310 Wilks)
       csihash['POD'] = matchra.sum() / (matchra.sum() + onlyra1.sum())
       # false alarm ratio (p 310 Wilks)
       csihash['FAR'] = onlyra2.sum() / (matchra.sum() + onlyra2.sum())
       csihash['PCORR'] = pcorr
       if verbose ==1:
          print('HERE' , matchra.size, matchra.shape , onlyra2.shape)
          print(matchra.sum() , onlyra1.sum() , onlyra2.sum())
          print('CSI POD FAR' , csihash['CSI'] , csihash['POD'] , csihash['FAR'])
 
    return csihash
