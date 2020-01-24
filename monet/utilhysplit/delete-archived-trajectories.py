#!/hysplit_users/allisonr/anaconda3/bin/python
#delete-archived-trajectories.py
#Deletes old trajectory files created from VOLCAT xml alert files
#from archive folder (older than 30 days)
#Generates log of activity: ~/VOLCANO/DeleteCron.log
import sys
import re
import os
import glob
import shutil
import logging
from datetime import datetime, timedelta

# start editable vars #
days_old = 20	# how old the files have to be before they are deleted
original_folder = "/hysplit-users/allisonr/AlertArchive/Trajectory/"
logfile = "/hysplit-users/allisonr/AlertArchive/DeleteTraj.log"	# recording deletion
# end editable vars #

# start function definitions #
def log(level,msg,tofile=True):
	print( msg)
	if tofile == True:
		if level == 0:
			logger.info(msg)
		else:
			logger.error(msg)
			
def end(code):
	log(0,"End.")
	log(0,"-------------------------")
	sys.exit(code)
# end function definitions #

# start process #
deletedate = datetime.now() - timedelta(days = days_old)
#Create list of files
allfiles = []
for name in glob.glob(original_folder+'*_[0-9]*'):
	allfiles.append(name)
#Create list of dates
dates = []
for x in allfiles:
	dates.append(re.search("(\d{14})", x).group(1))
#Create datetime object from extracted date string
dates2 = [datetime.strptime(x, '%Y%m%d%H%M%S') for x in dates]
#Create list of files to move
deletefiles = []
dd = 0
while dd < len(allfiles):
	if (dates2[dd] < deletedate):
		deletefiles.append(allfiles[dd])
	dd += 1

logger = logging.getLogger("cuarch")
hdlr = logging.FileHandler(logfile)
hdlr.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
logger.addHandler(hdlr) 
logger.setLevel(logging.INFO)

log(0,"Initializing...")

count = 0
size = 0.0
mm = 0
while mm < len(deletefiles):
    srcfile = deletefiles[mm]
    size = size + (os.path.getsize(srcfile) / (1024*1024))
    os.remove(srcfile)
    log(0,"Deleted '" + deletefiles[mm] + "'.")
    count = count + 1
    mm += 1

log(0,"Deleted " + str(count) + " files, totalling " + str(round(size,2)) + " MB.")
end(0)
