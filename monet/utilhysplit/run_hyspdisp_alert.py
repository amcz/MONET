#run_hyspdisp_alert.py
#reads volcat alert xml file, creates hysplit disp. control files, renames setup file
#runs hysplit with created control file
#makes concplot from cdump file

def make_files( fpath, fname):
    """ Makes CONTROL and SETUP files for Dispersion simulation from information in
    VOLCAT alert xml files pushed to the /pub/jpsss_upload/ ftp folder.
    
    Met Model set: GFS0p25
    Working Directory set: /hysplit-users/allisonr/VOLCANO/Dispersion/ - should be changed
    Variables set: duration (6 hours - forward)
    Vertical mixing scheme: use data - defaults to hysplit default

    These variables can be changed if desired. """

    from monet.util import volcalert
    from monet.utilhysplit import hcontrol
    from monet.utilhysplit import forecast_data as fd
    from monet.util import read_csv as usgs
    from datetime import timedelta as td
    from datetime import datetime as dtime
    import numpy as np
    import os

    met_model = 'GFS0p25' 
    wdir = '/hysplit-users/allisonr/VOLCANO/Dispersion/'
    duration = 6 #For dispersion run only
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
    pid = datetim+'_'+wmo_id

    #Read in list of Volcanoes and find vent height of closest volcano
    #used for calculating altitude
    csv = usgs.open_file(fname = '/hysplit-users/allisonr/Alice/MONET/monet/data/usgs_table.csv',delimiter = ',')
    headers = usgs.find_headers(csv)
    csv[headers[4]] = csv[headers[4]].str.normalize('NFKD').str.encode('ascii', errors = 'ignore').str.decode('utf-8')
    
    if csv[headers[4]].str.contains(volcname).any():
        volcname2 = volcname
        this_volc = csv[csv[headers[4]].str.contains(volcname)]
        alt = this_volc[headers[12]].values[0]
    else:
        volcname2 = 'Unknown'
        alt = 300.0
    if volcname.find('[') != -1:
        volcname2 = volcname[:volcname.index('[')]
        this_volc = csv[csv[headers[4]].str.contains(volcname2)]
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
        edate = sdate + td(days = ndays)
        #met = fd.findcycles_archive(sdate,edate,met_model,'Forward')
        met = fd.findcycles_forecast(sdate,met_model)
    if duration < 0:
        ndays = (duration // 24)
        edate = sdate + td(days = ndays)
        met = fd.findcycles_archive(sdate,edate,met_model,'Back')

    #Get more info from Alert
    hgts = volcalert.get_height(alerts)
    mass = volcalert.get_mass(alerts)
    radius = volcalert.get_radius(alerts)
    area = volcalert.get_area(alerts)

    #Write CONTROL file IF ASH ALERT
    control = hcontrol.HycsControl(fname='CONTROL.004d', working_directory = wdir)
    control.read()
    control.concgrids[0].outdir = wdir
    control.concgrids[0].outfile = 'cdump.'   #adds the pid automatically
    control.concgrids[0].sample_start = sdate.strftime("%y %m %d %H %M")
    levels = [alt, alt+3000, alt+10000]
    control.concgrids[0].levels = levels
    control.concgrids[0].nlev = np.size(levels)
    control.rename('CONTROL.disp' + str(pid), working_directory = wdir)
    for cgrid in control.concgrids:
        cgrid.outfile += str(pid)
    #simulation start date
    control.add_sdate(sdate)
    #set location of volcano
    control.remove_locations()
    control.add_location(latlon = (start_lat, start_lon), alt = alt)
    control.add_location(latlon = (start_lat, start_lon), alt = (float(hgts[0])*1000))
    #set the duration of the simulation (hours)
    control.run_duration = duration
    #set the duration of the emission
    for spec in control.species:
        spec.duration = 0.5
        # add met files
    control.remove_metfile(rall = True)
    for metdir, metfile in zip(met[0],met[1]):
        control.add_metfile(metdir, metfile)
    control.write(annotate = False)

    #Rename SETUP File for Dispersion
    setup = hcontrol.NameList(fname = 'SETUP.004d', working_directory = wdir)
    setup.read()
    setup.rename(name = 'SETUP.disp'+ str(pid), working_directory = wdir)
    setup.add('poutf', "'"'PARDUMP.disp'+str(pid)+"'")
    setup.write()

    #Writes a MAPTEXT file, used to add more information about each eruption to gifs.
    #Should have file nameing convention: MAPTEXT.'pid'
    #Contains lines of text
    with open(wdir + 'MAPTEXT.CFG', 'w') as fid:
        fid.write('Line1 \n')
        fid.write('Line2 \n')
        fid.write('HYSPLIT Dispersion Calculation for '+volcname2+' Volcano\n')
        fid.write('Line3 \n')
        fid.write('\n')
        fid.write('Line5 \n')
        fid.write('Line6 \n')
        fid.write('1-hr Average at the 6-hour Forecast of Volcanic Ash\n')
        fid.write('Lat: '+str(start_lat)+'    Lon: '+str(start_lon)+'    Vent Height: '+str(alt)+'m     Ash Column Height: '+str(float(hgts[0])*1000)+'m\n')
        fid.write('Ash Emission for 0.5 hours \n')
        fid.write('Alert Type: '+data[0]+'    VAAC: '+str(vaac)+'\n')
        fid.write('\n')
        current=dtime.utcnow()
        fid.write('Job Start: '+current.strftime('%B %d, %Y')+' at '+current.strftime('%H:%M:%S')+' UTC \n')
        fid.close()

    return wdir, pid, alert_type

def run_hysp(wdir, pid):
    """ Launches HYSPLIT using CONTROL and SETUP files created in make_files."""
    import os
    os.chdir(wdir)
    os.system('/hysplit-users/allisonr/HYSPLIT/exec/hycs_std  disp' +pid)
    return 'HYSPLIT Dispersion run for '+pid+' is completed!'

def make_disp(wdir, pid):
    """ Runs concplot to create figure from the cdump file created from run_hysp."""
    import os
    os.chdir(wdir)
    os.system('/hysplit-users/allisonr/HYSPLIT/exec/concplot -icdump.'+str(pid)+' -odisp_'+str(pid)+'.ps -a3 -c0 -d1 -e3 -f1 -r1 -s0 -11E-17 -51 -80 -91 +n')
    os.system('cp MAPTEXT.CFG MAPTEXT.disp'+pid)
    return 'Figure called disp_'+pid+'****.ps is created!'
    
def compress_kml(wdir, pid):
    import os
    """ Compresses kml files generated by concplot into kmz files. """
    os.chdir(wdir)
    os.system('zip disp_'+pid+'.kmz disp_'+pid+'_ps.kml')
    return 'kmz file called disp_'+pid+'.kmz is created!'

def check_for_figs(wdir, pid):
    import os
    """ Check to see if dispersion figures were created. """
    os.chdir(wdir)
    concfile = os.path.isfile('disp_'+pid+'0001.ps')
    concfile2 = os.path.isfile('disp_'+pid+'0002.ps')
    concfile3 = os.path.isfile('disp_'+pid+'0003.ps')
    return concfile, concfile2, concfile3

