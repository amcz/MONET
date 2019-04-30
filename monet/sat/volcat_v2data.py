#volcat_monet.py
#/hysplit-users/allisonr/VOLCANO/
#My attempt at creating a reader for VOLCAT data using xarray
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

#ModisHDF class -  gets satellite retrieval information from HDF file.

class VolcatHDF(object):
      """reads data from a hdf file"""
      def __init__(self, fname, verbose=0, datatype='532'):
            self.fname = fname
            self.missing = -999.0  
            self.version=[]
            self.lat=[]
            self.lon=[]
            self.height=[]
            self.radius=[]
            self.mass=[]
            self.dset=[]
            self.get_data(verbose=verbose)
            #self.rval=[]
            #self.mval=[]
            #self.get_latlon
            #self.hval=[]

      def get_data(self, verbose=1):
            """retrieves data from VOLCAT files."""
            self.dset = xr.open_dataset(self.fname,mask_and_scale=False,decode_times=False)
            self.version=self.dset.attrs['Default_Ash_Dust_Retrieval_Name']
            self.height = self.dset[self.version + '_ash_top_height'][:,:]
            self.radius = self.dset[self.version + '_ash_effective_radius'][:,:]
            self.mass = self.dset[self.version + '_ash_mass_loading'][:,:]
            lat = self.dset['pixel_latitude']
            lon = self.dset['pixel_longitude']
            self.lat = lat.attrs['scale_factor'] * lat[:,:]
            self.lon = lon.attrs['scale_factor'] * lon[:,:]
            if verbose:
               print('ATTRIBUTES')
               print(self.dset.attrs)
               print('VARIABLES')
               print(self.dset)

      def get_latlon(self):
            return  self.lat, self.lon
      
      def get_radius(self,unit='um'):
            """Returns 2d array of ash effective radius"""
            """Default units are micrometer"""
            rval = self.radius.where(self.radius != self.missing)
            return rval

      def get_mass_loading(self,unit='g/m2'):
            """Returns array with retrieved ash mass loading"""
            """Default units are grams/square meter"""
            mval = self.mass.where(self.mass != self.missing)
            return mval

      def get_height(self,unit='km'):
            """Returns array with retrieved height of the highest layer of ash."""
            """Default units are km above sea-level"""
            hval = self.height.where(self.height != self.missing)
            return hval

      def quickplot(self): 
          fig = plt.figure(1)
          val1=self.get_height()
          lat,lon = self.get_latlon()
          cb1 = plt.pcolormesh(lon,lat, val1)
          plt.colorbar(cb1)
          plt.title('Ash top height (km)')
          plt.xlabel('Longitude')
          plt.ylabel('Latitude')
          fig = plt.figure(2)
          val2=self.get_mass_loading()
          cb2 = plt.pcolormesh(lon,lat, val2)
          plt.colorbar(cb2)
          plt.title('Ash mass loading (g/m^2)')
          plt.xlabel('Longitude')
          plt.ylabel('Latitude')
          fig = plt.figure(3)
          val3=self.get_radius()
          cb3 = plt.pcolormesh(lon,lat, val3)
          plt.colorbar(cb3)
          plt.title('Ash effective radius (um)')
          plt.xlabel('Longitude')
          plt.ylabel('Latitude')
          plt.show()

      def geomap_plot(self):
            fig = plt.figure('Ash Top Height')
            val1=self.get_height()
            lat,lon = self.get_latlon()
            m=plt.axes(projection=ccrs.PlateCarree())
            m.add_feature(cfeat.LAND)
            m.add_feature(cfeat.COASTLINE)
            m.add_feature(cfeat.BORDERS)
            plt.pcolormesh(lon,lat, val1,transform=ccrs.PlateCarree())
            plt.colorbar()
            plt.title('Ash top height (km)')

            fig = plt.figure('Ash Mass Loading')
            val2=self.get_mass_loading()            
            m2=plt.axes(projection=ccrs.PlateCarree())
            m2.add_feature(cfeat.LAND)
            m2.add_feature(cfeat.COASTLINE)
            m2.add_feature(cfeat.BORDERS)
            plt.pcolormesh(lon,lat, val2,transform=ccrs.PlateCarree())
            plt.colorbar()
            plt.title('Ash mass loading (g/m^2)')

            fig = plt.figure('Ash Effective Radius')
            val3=self.get_radius()            
            m3=plt.axes(projection=ccrs.PlateCarree())
            m3.add_feature(cfeat.LAND)
            m3.add_feature(cfeat.COASTLINE)
            m3.add_feature(cfeat.BORDERS)
            plt.pcolormesh(lon,lat, val3,transform=ccrs.PlateCarree())
            plt.colorbar()
            plt.title('Ash effective radius (um)')

            plt.show()
            
