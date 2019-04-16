#from monet.util import volcatv2_data
import monet.util.tools
import monet.util.volcat_v3data as volcat
import matplotlib.pyplot as plt
import monet.models.hysplit as hysplit
import datetime
import seaborn as sns
import numpy as np
sns.set()


def sum_species(dset, plist, time, z, fillz=False):
    """
    dset : data array
    plist : list of names of species to add
    z   : level
    time : time
    fillz : boolean
    """
    returnra = []
    iii=0
    for species in plist:
        par = dset[species]
        val2d = (par.isel(time=time, z=z))
        if fillz : val2d = val2d.fillna(0)
        val2d = np.array(val2d)
        if iii==0:
            returnra = val2d
        else:
            returnra += val2d
    return returnra


class ConcPlot(object):
    
    def __init__(self, fname, d1=None, d2=None):
       self.fname = fname
       self.d1 = d1
       self.d2 = d2
       self.hxr = hysplit.open_dataset(fname, drange=[d1,d2])

       lats = self.hxr.coords['latitude']
       lons = np.array(self.hxr.coords['longitude'])
       self.latitude =  lats[0] 
       self.longitude = np.transpose(lons)[0]
       alist = self.hxr.data_vars.keys()
       print('HERE')
       print(list(alist))
       alist = list(alist)
       temp = self.hxr[alist[0]] 
      

    def getra(self, lev, tm, splist=None, fillz=False):
        """
        would like to make this so input time would be a datetime object
        and lev would be a height in meters.
 
        currently these are the indices.
        """
        if not splist: splist = list(self.hxr.data_vars(keys())
        val2d = sum_species(self.hxr, splist, time=tm,
                                    z=lev, fillz=fillz)
        return val2d


    def stepthru(self, splist=None, fillz=False):
        """
        splist: list of strings
        sums over list of species in the list.
        if None will sum over all species present in the file.
        """
        contour = False
        if not splist: splist = list(self.hxr.data_vars(keys())
        #hxr = hysplit.open_dataset(fname, drange=[d1,d2])
        #print('DATASET ************************')
        #print(hxr)
        #print('PAR1 ************************')
        #print(hxr['par1'])
        for lev in range(0,self.hxr.dims['z']): 
            for tm in range(0,self.hxr.dims['time']): 
                val2d = sum_species(self.hxr, splist, time=tm,
                                    z=lev, fillz=fillz)
                yield lev, tm, val2d

    def plotall(self, contour=False):
       longitude = self.longitude
       latitude = self.latitude
       splist = list(self.hxr.data_vars.keys())
       for lev, tm, val2d in self.stepthru(splist, fillz=contour):
           if not contour: cb = plt.pcolormesh(longitude, latitude, np.transpose(val2d))
           if contour : cb = plt.contourf(longitude, latitude, np.transpose(val2d))
           tstr = str(tm)
           plt.title('Time ' + tstr + ' Level ' + str(lev))
           plt.colorbar(cb)
           plt.show()
        


