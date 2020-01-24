#run_hysp_alert.py
#reads volcat alert xml file, creates hysplit traj. control files, renames setup file
#runs hysplit with created control file
#makes trajplot from tdump file
#IF tajectory plot figures are not created, the alert filename is not added to the list
#and the trajectory is initiated again at the next call of this program

def make_files( fpath, fname):
    """ Makes CONTROL and SETUP files for Trajectory simulation from information in
    VOLCAT alert xml files pushed to the /pub/jpsss_upload/ ftp folder.
    
    Met Model set: GFS0p25
    Working Directory set: /hysplit-users/allisonr/VOLCANO/Trajectory/ - should be changed
    Variables set: duration (36 hours - forward)
    Topbound set: top of model domain (m-agl)
    Vertical mixing scheme: use data - defaults to hysplit default

    These variables can be changed if desired. """

    from monet.util import volcalert
    from monet.utilhysplit import hcontrol
    from monet.utilhysplit import forecast_data as fd
    from monet.util import read_csv as usgs
    from datetime import timedelta as td
    from datetime import datetime as dtime
    import os

    met_model = 'GFS0p25' 
    wdir = '/hysplit-users/allisonr/VOLCANO/Trajectory/'
    duration = 36 #Forward: (+hours) Backward: (-hours)
    topbound = 30000  #Top of model domain (m-agl)
    vertical = 0 # 0:data, 1:isob, 2:isen, 3:dens 4:sigma 5:diverg 6:msl2agl 7:avg

    alerts = volcalert.open_xml(fpath+fname)
    #Returns an 8 dimensional array
    data = volcalert.get_alert_values(alerts)
    alert_type = data[1]
    start_lat = data[5]
    start_lon = data[6]
    #Returns the nearest volcanoes to the alert as list of 6 arrays
    nearby = volcalert.get_nearby(alerts)
    #Returns closest volcano 
    #Calculated based on minimum distance - requires all nearby values as input
    closest = volcalert.get_closest(nearby[0], nearby[1], nearby[2], nearby[3], nearby[4], nearby[5])
    volcname = closest[0]
    datetim = data[7].strftime('%Y%m%d%H%M%S')
    wmo_id = data[8]
    vaac = volcalert.get_vaac(alerts)
    pid = datetim+'_'+closest[4]

    #Read in list of Volcanoes and find vent height of closest volcano
    #used for calculating altitude
    csv = usgs.open_file(fname = '/hysplit-users/allisonr/Alice/MONET/monet/data/usgs_table.csv',delimiter = ',')
    headers = usgs.find_headers(csv)
    csv[headers[4]] = csv[headers[4]].str.normalize('NFKD').str.encode('ascii', errors = 'ignore').str.decode('utf-8')
    
    #If the volcano name has a "[" in the name (two names possible)
    if volcname.find('[') != -1:
        volcname2 = volcname[:volcname.index('[')]
        this_volc = csv[csv[headers[4]].str.contains(volcname2)]
        alt = this_volc[headers[12]].values[0]
    else:
        if csv[headers[4]].str.contains(volcname).any():
            volcname2 = volcname
            this_volc = csv[csv[headers[4]].str.contains(volcname)]
            alt = this_volc[headers[12]].values[0]
        else:
            volcname2 = 'Unknown'
            alt = 300.0

    #Find met fields necessary for the hysplit run
    #Two possible functions to gather appropriate met fields:
    #      findcycles_archive(datetime object start, datetime object end, 
    #                                    met model,'Forward' or 'Back')
    #      findcycles_forecast(datetime object start, met model)
    sdate = data[7]
    #Number of days for met fields needed - based on hour of trajectories
    if duration > 0:
        ndays = (duration // 24) + 1
        edate = sdate + td(days = ndays)
        #met = fd.findcycles_archive(sdate,edate,met_model,'Forward')
        met = fd.findcycles_forecast(sdate,met_model)
    if duration < 0:
        ndays = (duration // 24)
        edate = sdate + td(days = ndays)
        met = fd.findcycles_archive(sdate,edate,met_model,'Back')

    #Write CONTROL file
    control = hcontrol.HycsControl(fname='CONTROL.' + str(pid), working_directory = wdir, rtype = 'trajectory')
    # simulation start date
    control.add_sdate(sdate)
    # set the duration of the simulation (hours)
    control.run_duration = duration
    # set location of volcano
    control.add_location(latlon = (start_lat, start_lon), alt = alt, rate = False, area = False)
    control.add_location(latlon = (start_lat, start_lon), alt = alt+2000, rate = False, area = False)
    control.add_location(latlon = (start_lat, start_lon), alt = alt+4000, rate = False, area = False)
    control.add_location(latlon = (start_lat, start_lon), alt = alt+6000, rate = False, area = False)
    control.add_location(latlon = (start_lat, start_lon), alt = alt+10000, rate = False, area = False)
    control.add_location(latlon = (start_lat, start_lon), alt = alt+13000, rate = False, area = False)
    control.add_location(latlon = (start_lat, start_lon), alt = alt+16000, rate = False, area = False)
    control.vertical_motion = vertical
    control.ztop = topbound
    control.outdir = wdir
    control.outfile = 'tdump.'+str(pid)
    # add met files
    for metdir, metfile in zip(met[0],met[1]):
        control.add_metfile(metdir, metfile)
    control.write(annotate = False)

    #Rename SETUP file
    setup=hcontrol.NameList(fname='SETUP.Reventador.20190425', working_directory=wdir)
    setup.read()
    setup.rename(name='SETUP.' + str(pid), working_directory=wdir)
    setup.write(verbose=False)

    #Writes a MAPTEXT file, used to add more information about each eruption to gifs.
    #Should have file nameing convention: MAPTEXT.'pid'
    #Contains lines of text
    with open(wdir + 'MAPTEXT.CFG', 'w') as fid:
        fid.write('HYSPLIT Trajectory Calculation for '+volcname2+' Volcano\n')
        fid.write('Line2 \n')
        fid.write('Lat: '+str(start_lat)+'    Lon: '+str(start_lon)+'    Vent Height: '+str(alt)+' m\n')
        fid.write('Alert Type: '+data[0]+'    VAAC: '+str(vaac)+'\n')
        current=dtime.utcnow()
        fid.write('Job Start: '+current.strftime('%B %d, %Y')+' at '+current.strftime('%H:%M:%S')+' UTC \n')
        fid.close()

    return wdir, pid, alert_type

def run_hysp(wdir, pid):
    """ Launches HYSPLIT using CONTROL and SETUP files created in make_files.
    Prints current working directory."""
    import os
    cwd=os.getcwd()
    os.chdir(wdir)
    cwd2=os.getcwd()
#    print('Current wdir is: '+cwd2)
    os.system('/hysplit-users/allisonr/HYSPLIT/exec/hyts_std ' +pid)
    return 'HYSPLIT run for '+pid+' is completed!'

def make_traj(wdir, pid):
    """ Runs trajplot to create figure from the tdump file created from run_hysp."""
    import os
    
    cwd=os.getcwd()
    os.chdir(wdir)
    #os.system('python /hysplit-users/allisonr/Alice/MONET/monet/utilhysplit/hysplit_graf/trajplot.py -itdump.'+pid+' -otraj_'+pid+'.ps -a3 -A1 -l6 -s1 -v1 -k2 +n')
    os.system('/hysplit-users/allisonr/HYSPLIT/exec/trajplot -itdump.'+pid+' -otraj_'+pid+'.ps -a3 -A1 -l6 -s1 -v1 -k2 +n')
    os.system('cp MAPTEXT.CFG MAPTEXT.'+pid)
    return 'Figure called traj_'+pid+'.ps is created!'
    
def compress_kml(wdir, pid):
    import os
    """ Compresses kml files generated by trajplot into kmz files. """
    cwd = os.getcwd()
    os.chdir(wdir)
    os.system('zip traj_'+pid+'.kmz traj_'+pid+'_01.kml')
    return 'kmz file called traj_'+pid+'.kmz is created!'

def make_gif(wdir, pid):
    import os
    """ Converts ps files generated by trajplot into gif files. """
    cwd = os.getcwd()
    os.chdir(wdir)
    os.system('convert -antialias -flatten -trim -density 150 -resize 60% traj_'+pid+'.ps traj_'+pid+'.gif')
    return 'gif file called traj_'+pid+'.gif is created!'

def check_for_figs(wdir, pid):
    import os
    """ Check to see if trajectory figure was created. """
    cwd = os.getcwd()
    os.chdir(wdir)
    trajfile = os.path.isfile('traj_'+pid+'.gif')
    return trajfile

