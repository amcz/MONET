#volcat.py
#/hysplit-users/allisonr/VOLCANO/
#My attempt at creating a reader for VOLCAT data using xarray
#For use with MONET
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd

def open_dataset(fname):
      print(fname)
      dset = xr.open_dataset(fname,mask_and_scale=False,decode_times=False)
      #dset = _get_latlon(dset)
      #dset = _get_time(dset)
      #print(dset.latitude.scale_factor)
      return dset

def open_mfdataset(fname):
      #print(fname)
      #dset = xr.open_mfdataset(fname,concat_dim='time',decode_times=False,mask_and_scale=False)
      from glob import glob
      from numpy import sort
      files = sort(glob(fname))
      das = []
      for i in files:
            das.append(open_dataset(i))
      dset = xr.concat(das,dim='time') 
      #print(dset.latitude.scale_factor)
      #dset = _get_latlon(dset)
      return dset

def _get_time(dset):
      import pandas as pd
      test = str(dset.Image_Date)[1:] + str(dset.Image_Time)
      time = pd.to_datetime(test,format='%y%j%H%M%S')
      dset['time'] = time
      dset = dset.set_coords(['time'])
      dset = dset.expand_dims(dim='time')
      #dset.dims['time'] = time
      return dset

def _get_latlon(dset):
      lat = dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
      lon = dset.pixel_longitude[:,:] * dset.pixel_longitude.scale_factor
      dset = dset.rename({'pixel_latitude':'latitude'})
      dset = dset.rename({'pixel_longitude':'longitude'})
      dset = dset.set_coords(['latitude','longitude'])
      dset.latitude.load()
      dset.longitude.load()
      dset.coords['latitude'][:] = lat
      dset.coords['longitude'][:] = lon
      return dset
                  
def get_height(dset):      
      """Returns array with retrieved height of the highest layer of ash."""
      """Default units are km above sea-level"""
      var_header = dset.attrs['Default_Name_ash_ret']
      h_missing = dset[var_header + '_ash_top_height']._FillValue
      masked_height = ma.masked_equal(dset[var_header + '_ash_top_height'],h_missing)
      return masked_height

def get_radius(dset):
      """Returns 2d array of ash effective radius"""
      """Default units are micrometer"""
      var_header = dset.attrs['Default_Name_ash_ret']
      r_missing = dset[var_header + '_ash_effective_radius']._FillValue
      masked_radius = ma.masked_equal(dset[var_header + '_ash_effective_radius'],r_missing)
      return masked_radius

def get_mass_loading(dset):
      """Returns 2d array of ash effective radius"""
      """Default units are micrometer"""
      var_header = dset.attrs['Default_Name_ash_ret']
      m_missing = dset[var_header + '_ash_mass_loading']._FillValue
      masked_mass = ma.masked_equal(dset[var_header + '_ash_mass_loading'],m_missing)
      return masked_mass

def plot_height(dset):
       fig = plt.figure('Ash_Top_Height')
       #lat=dset.latitude
       #lon=dset.longitude
       lat = dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
       lon = dset.pixel_longitude[:,:]*dset.pixel_longitude.scale_factor
       masked_height=get_height(dset)
       m = plt.axes(projection=ccrs.PlateCarree())
       m.add_feature(cfeat.LAND)
       m.add_feature(cfeat.COASTLINE)
       m.add_feature(cfeat.BORDERS)
       plt.pcolormesh( lon, lat, masked_height, transform=ccrs.PlateCarree())
       plt.colorbar()
       plt.title('Ash Top Height (km)')
       plt.show()

def plot_radius(dset):
      fig = plt.figure('Ash_Effective_Radius')
      lat = dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
      lon = dset.pixel_longitude[:,:]*dset.pixel_longitude.scale_factor
      #lat=dset.latitude
      #lon=dset.longitude
      masked_radius=get_radius(dset)
      m=plt.axes(projection=ccrs.PlateCarree())
      m.add_feature(cfeat.LAND)
      m.add_feature(cfeat.COASTLINE)
      m.add_feature(cfeat.BORDERS)
      plt.pcolormesh( lon, lat, masked_radius, transform=ccrs.PlateCarree())
      plt.colorbar()
      plt.title('Ash effective radius (um)')
      plt.show()

def plot_mass(dset):
      fig = plt.figure('Ash_Mass_Loading')
      lat = dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
      lon = dset.pixel_longitude[:,:]*dset.pixel_longitude.scale_factor
      #lat=dset.latitude
      #lon=dset.longitude
      masked_mass=get_mass(dset)
      m=plt.axes(projection=ccrs.PlateCarree())
      m.add_feature(cfeat.LAND)
      m.add_feature(cfeat.COASTLINE)
      m.add_feature(cfeat.BORDERS)
      plt.pcolormesh( lon, lat, masked_mass, transform=ccrs.PlateCarree())
      plt.colorbar()
      plt.title('Ash mass loading (g/m^2)')
      pltshow()
