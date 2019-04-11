#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *
import numpy as np
import numpy.ma as ma
import matplotlib.dates as mdates
import string
import datetime
from scipy import interpolate
import subprocess
import sys, getopt
from operator import itemgetter
#import codecs
import os.path
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

##Collection of useful function


##removefromset

##autocorr - computes autocorrelation. 
  
##um2phi - converts grain size in microns to grain size in phi

##phi2um - converts phi to microns

##cat - concatenates files

##isnan - finds nans

##crop_contour_plot(x,y,z, xr=[], yr=[]):

##cropra(dra , missing, maskthera=0, crop=1, datara=[], buffer=0):

##histoplot 

def crop_contour_plot(x,y,z, xr=[], yr=[]):
    """crops 2d arrays so that x data lies in range given by xr and y
       data lies in range given by yr. xr and yr need to be list/array with [xmin, xmax]"""
    print 'YR', yr
    print 'XR', xr
    if yr != []:
        vp = np.where(np.logical_and(y> yr[0] , y < yr[1]))
        min_yi = np.min(vp[1])
        max_yi = np.max(vp[1]) 
        print min_yi , max_yi
        x=x[0:,min_yi:max_yi]
        y=y[0:,min_yi:max_yi]
        z=z[0:,min_yi:max_yi]
    if xr != []:
        vp = np.where(np.logical_and(x> xr[0] , x < xr[1]))
        min_xi = np.min(vp[0])
        max_xi = np.max(vp[0]) 
        print min_xi , max_xi
        x=x[min_xi:max_xi,]
        y=y[min_xi:max_xi,]
        z=z[min_xi:max_xi,]
    return x , y , z 



def cropra(dra , missing, maskthera=0, crop=1, datara=[], buffer=0):
    """Crops a 2d array to  rectangular area which contains all valid values."""
    """Set maskthearea=1 to return a masked array"""
    """if crop=0 then will mask but not crop array"""
    """if a datara is input then the datara will be cropped according to the valid area in xra"""
    if crop == 1:
	    vp = np.where(dra != missing)  #Indices of valid points
	    szx = dra.shape[0]
	    szy = dra.shape[1]
            ##if the min index < 0 Python seems to set it to zero
            ##if the max index > the size of the row/column then Python seems to set it to max.
            if vp[0] != [] and vp[1]!=[]:
		    if  np.min(vp[0]) - buffer >=0:
		        min_xi = np.min(vp[0]) - buffer
                    else:
                        min_xi = 0
		    max_xi = np.max(vp[0]) + buffer
		    if np.min(vp[1]) - buffer >=0:
		       min_yi = np.min(vp[1]) - buffer
                    else:
                       min_yi=0 
		    max_yi = np.max(vp[1]) + buffer
		    if datara == []:
		       cropdra = dra[min_xi:max_xi+1,min_yi:max_yi+1]
		    else:
		       cropdra = datara[min_xi:max_xi+1,min_yi:max_yi+1]
		    if maskthera == 1:
		       returnra = ma.masked_where(cropdra==missing , cropdra)
		    else:
		       returnra = cropdra
            else:   
                    if datara == []:
                       returnra = dra
                    else:
                       returnra = datara
                    print 'Data could  not be cropped. No valid data in array.'
    else: 
            returnra = ma.masked_where(dra== missing,dra)
    return returnra

def histoplot(x, bins=100, missing=(), pdf=1 , cdf=0, clr='b' , transparency=1, label='', 
              takelog=0, xrange=[]):
    xstats={}
    """creates 1d histogram. Returns dictionary with mean and median values."""
    align='mid'         #other values could be 'left' or 'right'
    histtype='bar'      #other values could be 'barstacked', 'step' or 'stepfilled'
    if missing != ():
       print 'MISSING' , missing
       xvalid = np.where(x!=missing)
       x = x[xvalid]
    if takelog == 1:
       x = np.log10(x)
    if xrange == []:
       plt.hist(x , bins=bins ,  normed=pdf , cumulative=cdf, color=clr , 
             alpha=transparency , label=label)
    else:
       plt.hist(x , bins=bins ,  normed=pdf , cumulative=cdf, color=clr , 
             alpha=transparency , label=label, range=xrange)
    #plt.show()
    xstats['mean'] = np.mean(x) 
    xstats['median'] = np.median(x) 
    return xstats


def roundtime(dtime):
    """rounds time to nearest hour. Rounds up for 30 minutes"""
    if dtime.minute >= 30:
       dt = datetime.timedelta(hours = 1)
    else:
       dt = datetime.timedelta(hours = 0)

    temptime = dtime.replace(minute=0) + dt
    return temptime   


def emap(quantity, lat , lon, range=[], clevs=[], pclrs=[], type='mesh', MISSING=0, title=''):
    #rcParams['figure.figsize']=8, 10
    #rcParams['axes.labelsize']=12
    #rcParams['xtick.labelsize']=12
    #rcParams['ytick.labelsize']=12
    #<F3>rcParams['font.family']='serif'
    pad = 1
    if range==[]:
       lonmin = np.min(lon) - pad 
       lonmax = np.max(lon) + pad
       latmin = np.min(lat) - pad
       latmax = np.max(lat) + pad
    else:
       lonmin = range[0]
       lonmax = range[1]
       latmin = range[2]
       latmax = range [3]
    if latmin < -90:
       latmin = -90
    if latmax > 90:
       latmax = 90
    if lonmin < -360:
       lonmin = -360
    if lonmax > 360:
       lonmax = 360
    if abs(lonmax-lonmin) < 10:
       minc = 2
    elif abs(lonmax-lonmin) < 40:
       minc = 5
    elif abs(lonmax-lonmin) < 80:
       minc = 10
    else:
       minc = 20

    lono = (lonmax - lonmin) / 2.0
    lato = (latmax - latmin) / 2.0

    #lono=10
    lato = 0
    print 'MIN' , np.min(quantity)
    print 'MAX' , np.max(quantity)
    if clevs==[]:
       diff = np.max(quantity) - np.min(quantity)
       ninc = diff/10
       clevs = range(floor(np.min(quantity)) , ceil(np.max(quantity)) , ninc)

    map = Basemap(projection='cyl', lon_0=lono, lat_0=lato,  llcrnrlon=lonmin, urcrnrlon=lonmax, llcrnrlat=latmin , urcrnrlat=latmax )
    
    map.drawlsmask(land_color="#E0E0E0", ocean_color='white', lakes=True)
    qma = ma.masked_equal(quantity, MISSING) 
    x, y = map(lon, lat)
    #cs= map.contourf(x, y, quantity, clevs)
    if type=='mesh':
        cm = plt.get_cmap('rainbow')
        cs= map.pcolormesh(x, y, qma, cmap=cm, vmin=clevs[0], vmax=clevs[-1])
    elif type=='fcontour':
        cs= map.contourf(x, y, qma, clevs)
    map.drawmeridians(np.arange(0,360,minc), labels=[0,0,0,1])
    map.drawparallels(np.arange(-90,90,minc), labels=[1,0,0,0])
    cbar = map.colorbar(cs, location='right', pad="5%")
    if title != '':
       plt.title(title)
    return map


def lastday(date):
    """given a date, returns date of last day in that month"""
    dtest = datetime.datetime(date.year, date.month+1, 1, 0)
    dtest = dtest - datetime.timedelta(days=1)
    return(dtest) 


def list2str(lst, sep=', ' ):
    """converts a list to a string"""
    lstr = sep.join(map(str,lst))
    return lstr

def logscale_fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a,b)

def dra2list(lat, lon, mv):
    """takes 2d arrays of matched lat, lon pairs and turns it into a list of tuples with missing
       values removed"""
    q = []
    for i in range(0,lat.shape[0]):
        for j in range(0, lat.shape[1]):
            q.append((lat[i][j], lon[i][j]))

    q = filter(lambda x: x[1] != mv and x[0] != mv, q)
    return zip(*q)

def removefromset(x,y,mv):
    """input x, y data and a value mv. Removes all data points for
    which either x or y is equal to the mv value"""
    xt = zip(x,y)
    xt = filter(lambda x: x[1] != mv and x[0] != mv, xt)
    #xt = filter(lambda x: x[1] <= mv and x[0] <= mv, xt)
    return zip(*xt)

def crop2dra(x,y,zs, xr=[], yr=[]):
    """crops 2d arrays so that x data lies in range given by xr and y
       data lies in range given by yr. xr and yr need to be list/array with [xmin, xmax]
       zs is a dictionary with 2d data arrays which correspond to the x, y positions."""
    print 'YR', yr
    print 'XR', xr

    if yr != []:
        vp = np.where(np.logical_and(y> yr[0] , y < yr[1]))
        min_yi = np.min(vp[1])
        max_yi = np.max(vp[1]) 
        print min_yi , max_yi
        x=x[0:,min_yi:max_yi]
        y=y[0:,min_yi:max_yi]
        for key in zs.keys:
            zs[key]=zs[key][0:,min_yi:max_yi]
    if xr != []:
        vp = np.where(np.logical_and(x> xr[0] , x < xr[1]))
        min_xi = np.min(vp[0])
        max_xi = np.max(vp[0]) 
        print min_xi , max_xi
        x=x[min_xi:max_xi,]
        y=y[min_xi:max_xi,]
        for key in zs.keys:
            zs[key]=zs[key][min_xi:max_xi,]
    return x , y , zs

def autocorr(x, t, n=10, mv=''):
    """Assumes that times, t, are evenly spaced.
       mv is the value that indicates data is missing."""
    rt = []
    tlist = []
    nlist = []
    for lag in range(n):
        #lag = 1 
        print lag
        dt = t[lag] - t[0]
        dt = dt.seconds/60.0/60.0 + dt.days*24.0  #convert to hours    
        print t[lag] , t[0]  , dt
        if lag != 0:
           x1 = (x[0:-1*lag])    
        else:
           x1 = (x[0:])    
        x2 = (x[lag:])    
 
        ##get rid of missing values.
        xt = zip(x1,x2)
        xt = filter(lambda x: x[1] != mv and x[0] != mv, xt)
        npts = len(xt)
        x1 , x2 = zip(*xt)
        #################################################

        x1 = np.array(x1)
        x2 = np.array(x2)

        x1bar = np.mean(x1)
        x2bar = np.mean(x2)
        a = (x1-x1bar) * (x2-x2bar)
        b = (x1-x1bar)**2  
        c = (x2-x2bar)**2  
        r = a.sum() / np.sqrt(b.sum() * c.sum())

        rt.append(r)
        tlist.append(dt)
        nlist.append(npts)
    return rt, tlist , nlist


def um2phi(um):
    """converts grain size in phi to microns"""
    """phi = log_2(d/do) where do-1mm"""
    do = 1000.0 #reference size in microns.
    phi = -1 * log(um/do, 2)
    return(phi) 


def phi2um(phi):
    """converts grain size in phi to microns"""
    """phi = log_2(d/do) where do-1mm"""
    do = 1000.0 #reference size in microns.
    d = do * 2**(-1 * phi)
    return(d) 


def cat(infiles, outfile):
    """Input is two text files. Adds the infile to the end of the outfile"""
    success=1
    with open(outfile, "w") as myoutfile:
        my_cmd = ["cat"] + infiles
        try:
           subprocess.call(my_cmd, stdout=myoutfile)
        except:
           success=-1
           print "python cat error in mytools.py. Failed to concatenate files"
    return success

def nan_finder(y):

    return np.isnan(y), lambda z: z.nonzero()[0]
         
def set_date_axis(ax, dates):
    """picks correct tick marks to be displayed on the x axis
    """
    dt = dates[-1]-dates[0]
    print dt
    if dt.days >=3 and dt.days <=5:
       print dt
       ml=mdates.DayLocator()
       mml=mdates.HourLocator()
       daysFmt=mdates.DateFormatter('%d %H')
       dstr = 'Time (DD HH UTC in '
    elif dt.days > 5 and dt.days < 31:
       ml=mdates.DayLocator(bymonthday=range(0,31,5))
       mml=mdates.DayLocator()
       daysFmt=mdates.DateFormatter('%d')
       dstr = 'Time (Day in '
    elif dt.days > 31:
       ml=mdates.MonthLocator()
       mml=mdates.DayLocator(bymonthday=range(0,31,5))
       daysFmt=mdates.DateFormatter('%M %y')
       dstr = 'Time (Month Year) '
    elif dt.days < 3 and dt.days >=1:
       ml=mdates.HourLocator(byhour=range(0,24,6))
       mml=mdates.HourLocator()
       daysFmt=mdates.DateFormatter('%d %H')
       dstr = 'Time (DD HH UTC in '
    elif dt.seconds <= 6*60*60 :
       ml=mdates.HourLocator(byhour=range(0,24,2))
       mml=mdates.HourLocator()
       daysFmt=mdates.DateFormatter('%d %H')
       dstr = 'Time (DD HH UTC in '
    else:
       ml=mdates.HourLocator(byhour=range(0,24,2))
       mml=mdates.HourLocator()
       daysFmt=mdates.DateFormatter('%d %H')
       dstr = 'Time (DD HH UTC in '
    ax.xaxis.set_major_locator(ml)
    ax.xaxis.set_major_formatter(daysFmt)
    ax.xaxis.set_minor_locator(mml)
    year=dates[0].year
    month=dates[0].strftime("%b")
    ax.set_xlabel(dstr  + str(month) + '  ' + str(year) +')')

def sortxy(x, y):
   """Input x and y arrays which consitute a data set.
       Sorts by x value and returns re-ordered arrays."""
   xy = zip(x,y)
   xy = list(set(xy))   #remove duplicate xy pairs
   xy = sorted(xy, key=itemgetter(0))
   xy = zip(*xy)
   x = xy[0]
   y = xy[1]
   return x , y

def sortxyz(x, y, z):
   """Input x and y  and z arrays which consitute a data set.
       Sorts by z value and returns re-ordered arrays."""
   xyz = zip(x,y,z)
   xyz = sorted(xyz, key=itemgetter(2))
   xyz = zip(*xyz)
   x = xyz[0]
   y = xyz[1]
   z = xyz[2]
   return x , y , z

def spline_v(xra , yra , k=1):
    """return spline for each array input"""
    i=0
    outra=[]
    xra=np.atleast_2d(xra)
    yra=np.atleast_2d(yra)
    for x in xra:
        y=yra[i]
        #print x
        #print y
        f=interpolate.splrep(x,y,k=k) 
        outra.append(f)
        i+=1
    return  outra

    




def H2P(H):
    """convert height to pressure. If H in meters, returns P in pascals"""
    B=2.25577e-5
    A=101325.0
    C=5.25588
    P=A*(1-H*B)**C 
    return P 

def P2H(P):
    """convert pressure to height"""
    #P=np.array(pressure_list)
    B=2.25577e-5
    A=101325.0
    C=5.25588
    P=P*100   #convert from mbar to Pa
    h= 1/B * (1-(P/A)**(1/C))  
    #h=h.tolist()
    return h

def sphu2relh(sphu, p , t):
    """convert specific humidity to relative humidity"""
    R=287  #J/kg.K
    rho_air = p/(R*C2K(t))*1000 #convert from kg/m3 to g/m3
    VD = sphu * rho_air
    VDS = sat_vap_density(t)
    relh=sphu
    return(relh)

def sat_vap_density(t, type=1):
    """Calculates the saturated vapor density given a temperature"""
    if type==1:
        ##Find saturated vapor density at the temperature.
        ##use empirical fit for values of temp(oC) from 0 to 40.
        ##hyperphysics website http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/relhum.html#c3
        if t>=0:
            A=5.018
            B=0.32321
            C=8.1847e-3
            D=3.1243e-4
            VDS=A + B * t + C * t**2 + D*t**3  ##result in g/m3
        else: 
            VDS=0
    return VDS

def relh2sphu(relh, h , t):
    """convert relative humidity to specific humidity
     Not working correctly for higher altitutdes / temperatures below 0C.."""
    ##Find saturated vapor density at the temperature.
    ##use empirical fit for values of temp(oC) from 0 to 40.
    ##hyperphysics website http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/relhum.html#c3
    VDS= sat_vap_density(t)
    #VDS=saturated vapor density
    VD = VDS * relh / 100 ##actual vapor density g/m3
    ##Now need to find how how many grams of dry air are in a m3. 
    R=287  #J/kg.K
    rho_air = H2P(h)/(R*C2K(t))*1000 #convert from kg/m3 to g/m3
    if rho_air > 1300:
       rho_air = 1300
    rho_air = 1300 #convert from kg/m3 to g/m3
    sphu=VD/rho_air
    print  rho_air, VD , VDS , relh , t
    if sphu < 0 :
       print "warnging: negative sphu" , rho_air, VD , VDS , relh , t
       sphu=0
    return(sphu)

def is_number(s):
    try:
       float(s)
       return True
    except:
       return False


def addticks(ax, newlocs, newLabels, pos='x'):
    plt.draw()

    if pos=='x':
       locs = ax.get_xticks().tolist()
       labels = [x.get_text() for x in ax.get_xticklabels()]
        
    Dticks = dict(zip(locs, labels))
    
    for Loc, Lab in zip(newLocs,newLabels):
        Dticks[Loc] = Lab

    locs = list(Dticks.keys())
    labels = list(Dticks.values())

    if pos=='x':
       ax.set_xticks(locs)
       ax.set_yticks(labels)

     
 







