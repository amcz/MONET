from monet.utilhysplit.hcontrol import HycsControl 
from monet.utilhysplit.hcontrol import NameList
import datetime
import os


def make_setup(pid,
               wdir):
    if not os.path.isfile(wdir + 'SETUP.' + str(pid)):
        setup=NameList(fname='SETUP.VAAC', working_directory='./')
        setup.read()
        setup.rename(name='SETUP.' + str(pid), working_directory=wdir)
        # rename pardump file
        setup.nlist['poutf'] = setup.nlist['poutf'] + str(pid)
        #setup.write(verbose=False)
        return True
    else:
        return False

def make_control(pid,     
                 wdir, 
                 simulation_start,
                 runtime,
                 duration,
                 lat,
                 lon,
                 vent_ht,
                 plume_ht,
                 mdir,
                 mfiles
                 ):
    """
    Reads an existing CONTROL file and then modifies it and writes it out
    under a different name.
    pid : string (process id)
    wdir : string : directory where control file is to be written.
    simulation_start : datetime object
    runtime : duration of simulation
    duration : duration of emission
    lat, lon : start lat and lon
    vent_ht : vent height
    plume_ht : plume height
    mdir : met directories
    mfiles : met files
    Return
    gs : float concentration grid size spacing
    """
    control = HycsControl(fname='CONTROL.VAAC')
    control.read() 
    gs = control.concgrids[0].latdiff
    ftest=False
    if not os.path.isfile(wdir + 'CONTROL.' + str(pid)):
        ftest=True 
        # rename the control file
        control.rename('CONTROL.' + str(pid), working_directory=wdir)
        # add pid to all cdump outputs
        for cgrid in control.concgrids:
            cgrid.outfile +=  str(pid)
        # simulation start date
        control.add_sdate(simulation_start)
        # set the duration of the simulation (hours)
        control.run_duration = runtime
        # set the duration of the emission.
        for spec in control.species:
            spec.duration = duration 
        # set location of volcano
        control.remove_locations()
        control.add_location(latlon=(lat, lon), alt=vent_ht) 
        control.add_location(latlon=(lat, lon), alt=plume_ht) 
        # add met files
        control.remove_metfile(rall=True)
        for metdir, metfile in zip(mdir,mfiles):
            control.add_metfile(metdir, metfile)
        control.write(annotate=False) 
    return ftest, gs


def make_tcontrol(cname):
    """make control file for trajectory"""
    simulation_start = datetime.datetime.now()
    metdir='dir1'
    metfiles = 'file1'
    control = HycsControl(fname=cname, rtype='trajectory')
    control.add_location(latlon=(35, 50), alt=100) 
    control.add_location(latlon=(35, 40), alt=200) 
    control.add_sdate(simulation_start)
    control.add_duration=24
    control.outfile='tdump.01'
    control.add_metfile(metdir, metfiles)
    control.write()
#make_tcontrol('TEST')  

def make_traj_control(pid, wdir, simulation_start, runtime, lat, lon, height, mdir, mfiles):
    """ pid : string (process id)
         wdir : string : directory where control file is to be written.
         simulation_start : datetime object
         runtime : duration of simulation
         lat, lon : start lat and lon
         height : height(s) of trajectory start
         mdir : met directories
         mfiles : met files"""
    control = HycsControl(fname=wdir + 'CONTROL.' + str(pid),rtype = 'trajectory')
    # simulation start date
    control.add_sdate(simulation_start)
    # set the duration of the simulation (hours)
    control.run_duration = runtime
    # set location of volcano
    control.add_location(latlon=(lat, lon), alt=vent_ht) 
    control.add_location(latlon=(lat, lon), alt=plume_ht) 
    control.outfile = 'tdump.'+str(pid)
    # add met files
    for metdir, metfile in zip(mdir,mfiles):
        control.add_metfile(metdir, metfile)
    control.write()

