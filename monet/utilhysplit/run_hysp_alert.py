#run_hysp_alert.py
#reads volcat alert xml file, creates hysplit traj. control files, renames setup file
#runs hysplit with created control file
from monet.util import volcalert
from monet.utilhysplit import hcontrol
from monet.utilhysplit import forecast_data as fd
from monet.util import read_csv as usgs
from datetime import timedelta as td
import subprocess
import os

### VARIABLES NEED TO BE SET ###
#Path and name of VOLCAT alert xml file
fpath = '/hysplit-users/allisonr/VOLCANO/'
fname = 'geocatASHOBJ.GOES-16.Full.xml'
#Met model for trajectory
met_model = 'GDAS1'
#HYSPLIT working directory
wdir = '/hysplit-users/allisonr/VOLCANO/Trajectory/'
#Setting duration of simulation (hours), top of model domain (m), vertical motion method
duration = 36 #Forward: (+hours) Backward: (-hours)
topbound = 50000  #Top of model domain (m-agl)
vertical = 0 #0:data, 1:isob, 2:isen, 3:dens 4:sigma 5:diverg 6:msl2agl 7:avg
########################

#Opens the xml file
alerts = volcalert.open_xml(fpath+fname)
#Returns the alert_header, alert_type, status, confidence, method, lat of radiative center
#lon of radiative center, and datetime object of eruption start
# An 8 dimensional array
data = volcalert.get_alert_values(alerts)
start_lat = data[5]
start_lon = data[6]
#Returns the nearest volcanoes to the alert as separate  arrays of:
#Volcano name, latitude, longitude, distance from alert, thermal signature, and volc. ID
#6 arrays, listed 
nearby = volcalert.get_nearby(alerts)
#Returns the name, lat, lon, therm signature, and volc. ID of closest volcano 
#to the detected radiative signature returned from get_alert_values
#Calculated based on minimum distance - requires all nearby values as input
closest = volcalert.get_closest(nearby[0],nearby[1],nearby[2],nearby[3],nearby[4],nearby[5])
volcname = closest[0]
date = data[7].strftime('%Y%m%d')
pid = volcname + date

#Read in list of Volcanoes and find vent height of closest volcano - used for calculating altitude
csv = usgs.open_file(fname = '/hysplit-users/allisonr/Alice/MONET/monet/data/usgs_table.csv',delimiter = ',')
headers = usgs.find_headers(csv)
this_volc = csv[csv[headers[4]] == volcname]
alt = this_volc[headers[12]].values[0]

#Find met fields necessary for the hysplit run
#Two possible functions to gather appropriate met fields:
#      findcycles_archive(datetime object start, datetime object end, 
#                                    met model,'Forward' or 'Back')
#      findcycles_forecast(datetime object start, met model)
sdate = data[7]
#Number of days for met fields needed - based on hour of trajectories
if duration > 0:
    ndays = (duration // 24) + 1
    edate = data[7] + td(days = ndays)
    met = fd.findcycles_archive(sdate,edate,met_model,'Forward')
if duration < 0:
    ndays = (duration // 24)
    edate = data[7] + td(days = ndays)
    met = fd.findcycles_archive(sdate,edate,met_model,'Back')

#Write CONTROL file
control = hcontrol.HycsControl(fname='CONTROL.' + str(pid), working_directory = wdir, rtype = 'trajectory')
# simulation start date
control.add_sdate(sdate)
# set the duration of the simulation (hours)
control.run_duration = duration
# set location of volcano
control.add_location(latlon=(start_lat, start_lon), alt = alt, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+1000, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+1500, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+2000, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+3000, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+5000, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+7500, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+10000, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+12500, rate = False, area = False)
control.add_location(latlon=(start_lat, start_lon), alt = alt+15000, rate = False, area = False)
control.vertical_motion = vertical
control.ztop = topbound
control.outdir = wdir
control.outfile = 'tdump.'+str(pid)
# add met files
for metdir, metfile in zip(met[0],met[1]):
    control.add_metfile(metdir, metfile)
control.write(annotate = True)

#Rename SETUP file
setup=hcontrol.NameList(fname='SETUP.Reventador.20190425', working_directory=wdir)
setup.read()
setup.rename(name='SETUP.' + str(pid), working_directory=wdir)
setup.write(verbose=True)

#Run hysplit
cwd=os.getcwd()
print('Current wdir is: '+cwd)
os.chdir(wdir)
cwd2=os.getcwd()
print('Current wdir is: '+cwd2)
os.system('ls')
os.system('/hysplit-users/allisonr/hysplit/trunk/exec/hyts_std ' +pid)
