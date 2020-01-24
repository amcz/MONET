#!/hysplit_users/allisonr/anaconda3/bin/python
#move-old-files.py
#Moves old VOLCAT xml alert files out of ftp folder into archive folder
#Generates log of activity: MoveLog.log
import sys
import os
import glob
import shutil
import time
import logging
from datetime import datetime, timedelta

# start editable vars #
hours_old = 48	# how old the files have to be before they are moved
original_folder = "/pub/jpsss_upload/"	# folder to move files from
new_folder = "/hysplit-users/allisonr/AlertArchive/VOLCAT/"	# folder to move files to
logfile = "/hysplit-users/allisonr/VOLCANO/MoveLog.log"	# log file to record what has happened
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
move_date = time.mktime(movedate.timetuple())
logger = logging.getLogger("cuarch")
hdlr = logging.FileHandler(logfile)
hdlr.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
logger.addHandler(hdlr) 
logger.setLevel(logging.INFO)

log(0,"Initializing...")

count = 0
size = 0.0

for filename in glob.glob1(original_folder, "*.xml"):
    srcfile = os.path.join(original_folder, filename)
    destfile = os.path.join(new_folder, filename)
    if os.stat(srcfile).st_mtime < move_date:
        if not os.path.isfile(destfile):
            size = size + (os.path.getsize(srcfile) / 1024)
            shutil.move(srcfile, destfile)
            log(0,"Archived '" + filename + "'.")
            count = count + 1

log(0,"Archived " + str(count) + " files, totalling " + str(round(size,2)) + " KB.")
end(0)
