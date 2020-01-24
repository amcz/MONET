#!/hysplit_users/allisonr/anaconda3/bin/python
#move-old-trajectories.py
#Moves old trajectory files created from  VOLCAT xml alert files
#from trajectory output folder into archive folder
#Generates log of activity: MoveTrajs.log
import sys
import re
import os
import glob
import shutil
import logging
from datetime import datetime, timedelta

# start editable vars #
hours_old = 48	# how old the files have to be before they are moved
original_folder = "/hysplit-users/allisonr/VOLCANO/Trajectory/"	# from folder
new_folder = "/hysplit-users/allisonr/AlertArchive/Trajectory/"	# to folder
logfile = "/hysplit-users/allisonr/VOLCANO/MoveTraj.log"	# recording move
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
movedate = datetime.now() - timedelta(hours=hours_old)
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
movefiles = []
dd = 0
while dd < len(allfiles):
	if (dates2[dd] < movedate):
		movefiles.append(allfiles[dd])
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
while mm < len(movefiles):
    srcfile = movefiles[mm]
    size = size + (os.path.getsize(srcfile) / (1024*1024))
    shutil.move(srcfile, new_folder)
    log(0,"Archived '" + movefiles[mm] + "'.")
    count = count + 1
    mm += 1

log(0,"Archived " + str(count) + " files, totalling " + str(round(size,2)) + " MB.")
end(0)
