import xarray as xr
#import matplotlib.pyplot as plt
#import textwrap
#import cartopy.crs as ccrs
#import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
#from glob import glob
#from numpy import sort
#import pandas as pd
#import volcat
#import tools
from monet.models import hysplit
#import monet.util.csi_updates  as csi
#from monet.utilhysplit import hysp_func
#from monet.util import volcMER
import datetime
import seaborn as sns
import os
import monet.utilhysplit.hysp_func as hf
from netCDF4 import Dataset

# 01/28/2020 AMC cdump2awips created to make a netcdf file appropriate for input into AWIPS
# hysplit.py was modified in the forked version of MONET to make this work.

def combine_cdump(blist, d1=None,d2=None):
    # blist is a dictionary.
    # key is a tag indicating the source
    # value is tuple (filename, metdata)
    # this probably isn't the best way to arrange it.
    # possibly make a class to store that data.

    # d1 datetime object. first date to keep in DatArrayarray
    # d2 datetime object. last date to keep in DataArray

    # RETURNS
    # newhxr is an xarray data-array with 6 dimensions.
    #        lat, lon, time, level, ensemble tag, source tag
    # dt is the averaging time of the hysplit output.

    iii=0
    ylist = []
    dtlist = []
    sourcelist = [] 

    # first loop go through to get expanded dataset.
    xlist = []
    sourcelist = []
    enslist = []    
    for key in blist:
        xsublist = []
        for fname in blist[key]:
            if d1 and d2:
                hxr = hysplit.open_dataset(fname[0], drange=[d1,d2])
            else: #use all dates
                hxr = hysplit.open_dataset(fname[0])
            xrash = hf._add_species(hxr)
            xsublist.append(xrash)
            enslist.append(fname[1])
            dtlist.append(hxr.attrs['Sample time hours'])
            if iii==0:
               xnew = xrash.copy()
            else:
               a,xnew = xr.align(xrash, xnew)
            iii+=1 
        sourcelist.append(key)
        xlist.append(xsublist) 
    # xnew is now encompasses the area of all the data-arrays
    # now go through and expand each one to the size of xnew.
    iii=0
    jjj=0
    ylist = []
    slist = []
    for sublist in xlist:  
        hlist = []
        for temp in sublist:  
            # expand to same region as xnew
            aaa,bbb = xr.align(temp, xnew)
            aaa.expand_dims('ens')
            aaa['ens'] = enslist[iii] 
            iii+=1
            hlist.append(aaa)
        # concat along the 'ens' axis
        new = xr.concat(hlist, 'ens')
        print('HERE NEW', new)
        ylist.append(new)
        slist.append(sourcelist[jjj])
        jjj+=1
     
    dtlist = list(set(dtlist))
    print('DT', dtlist, dtlist[0])
    dt = dtlist[0]
    newhxr = xr.concat(ylist, 'source')
    print('sourcelist', slist)
    newhxr['source'] = slist
    print('NEW NEW NEW', newhxr)
    #newhxr['ens'] = metlist       

    # newhxr is an xarray data-array with 6 dimensions.
    # dt is the averaging time of the hysplit output.
    return newhxr, dt
        
def cdump2awips(flist, outname, format='NETCDF4', d1=None, d2=None):
    
    #hxr = hysplit.open_dataset(fname, drange=[d1,d2])
    xrash, dt = combine_cdump(flist, d1=d1, d2=d2)
    nra = xrash.values
    # array with top height of levels in the file.
    levelra = xrash.z.values
  
    sample_time = np.timedelta64(int(dt),'h')

    fid = Dataset(outname, 'w', format='NETCDF4')
    
    # GLOBAL ATTRIBUTES
    fid.SourceFiles = 'Kasatochi'
 
    print(xrash.shape)
    print(xrash.coords) 
    # DEFINE DIMENSIONS
    lat_shape= xrash.shape[4]
    lon_shape= xrash.shape[5]
    lat = fid.createDimension('latitude',lat_shape)
    lon = fid.createDimension('longitude',lon_shape)
    level = fid.createDimension('levels',len(levelra))

    # differnt runs with differnt sources
    source_shape=xrash.coords['source'].shape[0]
    source = fid.createDimension('source',source_shape)
    # differnt runs with different met data.
    ens_shape=xrash.coords['ens'].shape[0]
    ensemble = fid.createDimension('ensemble',ens_shape)

    time = fid.createDimension('time',1) # one time per file
    bnds = fid.createDimension('bnds',2) # two bounds per time.

    #origin = fid.createDimension('origins',hxr.attrs['Number Start Locations'])
    origin = fid.createDimension('origins',1)
    # Scalar variables
     
    #latra, lonra = hf.getlatlon(hxr)  
    latra = xrash.latitude[:,0]
    lonra = xrash.longitude[0]
    print('lat shape', latra.shape, xrash.shape)
    print('lon shape', lonra.shape, xrash.shape)
    # Define 2D variables with attributes
    concid = fid.createVariable('conc', 'f4', 
                               ('source', 'ensemble','levels','latitude','longitude'))
    print('SHAPE', source_shape, ens_shape, len(levelra), lat_shape, lon_shape) 
    concid.units = 'unit/m3'
    concid.long_name = 'Concentration Array - ash'
    #concid.meteorological_model_id = hxr['Meteorological Model ID']
    #concid.ensemble_member_tag = 1
    #concid.source_tag = 1

    # Dimension with different ensemble members.
    ensembleid = fid.createVariable('ensemble','str',('ensemble'))
    sourceid = fid.createVariable('source','str',('source'))

    latid = fid.createVariable('latitude', 'f4', ('latitude'))
    latid.long_name = 'latitude degrees north from the equator'
    latid.units = 'degrees_north'
    latid.point_spacing = 'even'
    lonid = fid.createVariable('longitude', 'f4', ('longitude'))
    lonid.long_name = 'longitude degrees east from the greenwhich meridian'
    lonid.units = 'degrees_east'
    lonid.point_spacing = 'even'

    levelid = fid.createVariable('levels','int',('levels'))
    # attributes for levels
    levelid.long_name = 'Top height of each layer'
    levelid.units='m'

    timeid = fid.createVariable('time', 'f4', ('time'))
    # attributes for time grid.
    timeid.units = 'days since 1970-01-01 00:00:00' 
    timeid.standard_name = 'time' 
    timeid.bounds = 'time_bnds' 
    timeid.calendar = 'gregorian'

    #olon = fid.createVariable('olon','f4',('origins'))
    #olat = fid.createVariable('olat','f4',('origins'))
    #olat.long_name = 'Simulation origin latitude'
    #olon.long_name = 'Simulation origin longitude'

    #olvl = fid.createVariable('olvl','f4',('origins'))
    #otim = fid.createVariable('otim','f4',('origins'))

    time_bnds = fid.createVariable('time_bnds','f4',('time','bnds'))

    # Put data into variables
    # only one time per file.
    epoch = np.datetime64('1970-01-01T00:00:00Z')
    date1 = xrash.time[0].values
    t1 = (xrash.time[0].values - epoch) / np.timedelta64(1,'s')
    # change seconds to days
    t1 = t1 / (24.0 * 60 * 60)
    
    t2 = ((xrash.time[0].values + sample_time) - epoch) / np.timedelta64(1,'s')
    t2 = t2 / (24.0 * 60 * 60)
  
    #oloc = hxr.attrs['Starting Locations']
      
    #oloc2 = list(zip(*oloc))
    #olat[:] = oloc2[0]
    #olon[:] = oloc2[1]
    temp = xrash.loc[dict(time=date1)]
    print(temp.values.shape) 
    print('date', date1)
    #concid[:] = xrash.loc[:,date1].values
    concid[:] = xrash.loc[dict(time=date1)].values
    latid[:] =  latra
    lonid[:] =  lonra
    levelid[:] = levelra
    timeid[:] = t1
    time_bnds[:] = [[t1 , t2]]

    ensembleid[:] = xrash.coords['ens'].values
    sourceid[:] = xrash.coords['source'].values
    fid.close()


def maketestblist():
    d1 = datetime.datetime(2008,8,8,12)
    d2 = datetime.datetime(2008,8,8,13)
    blist={}
    flist=[]
    dname='/pub/Scratch/alicec/KASATOCHI/cylindrical/e3/'
    fname='wrf.e3.bin'
    flist.append((os.path.join(dname,fname), 'WRF'))
    fname='gdas.e3.bin'
    flist.append((os.path.join(dname,fname), 'GDAS'))
    blist['S3'] = flist

    flist1=[]
    dname='/pub/Scratch/alicec/KASATOCHI/cylindrical/e2/'
    fname='wrf.e2.bin'
    flist1.append((os.path.join(dname,fname),'WRF'))
    fname='gdas.e2.bin'
    flist1.append((os.path.join(dname,fname),'GDAS'))
    blist['S2'] = flist1
    return blist

def maketestncfile():
    blist = maketestblist()
    oname = 'out.nc'
    d1 = datetime.datetime(2008,8,8,12)
    d2 = datetime.datetime(2008,8,8,13)
    cdump2awips(blist, oname, d1=d1, d2=d2)

def maketestra():
    #d1 = datetime.datetime(2008,8,8,10)
    #d2 = datetime.datetime(2008,8,8,13)
    d1=None
    d2=None
    blist = maketestblist()
    xrash, dt = combine_cdump(blist, d1=d1, d2=d2) 
    return xrash, dt

