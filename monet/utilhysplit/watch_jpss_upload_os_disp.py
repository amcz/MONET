#!/hysplit_users/allisonr/anaconda3/bin/python
#watch_jpss_upload_os.py
#Monitors ftp folder (/pub/jpsss_upload/) for new xml file
#If new file, then run file (run_hysp_alert.py)
#New file is found in directory: /hysplit-users/allisonr/Alice/MONET/monet/utilhysplit/
#Should create a hysplit control file, and run hysp trajectory whenever a new xml file 
#is detected in the ftp folder

import time
import os
import re
import json
import numpy as np
from monet.utilhysplit import run_hyspdisp_alert as rha2

#Path of ftp folder
path_to_watch = '/pub/jpsss_upload/'
#Creating list of current files in ftp folder
current = list(fi for fi in os.listdir(path_to_watch) if fi.endswith('.xml'))
#Reading the list of files from the last iteration
direct = '/hysplit-users/allisonr/VOLCANO/'
with open(direct+'original_list_disp.txt', 'r') as fs:
    original = json.loads(fs.read())
#Determining what was added and what was removed
added = [fs for fs in current if not fs in original]
removed = [fs for fs in original if not fs in current]

original1 = original
#If files were added, loop for all files added 
if added:
    print('\n')
    print('Added:', ',\n '.join(added),'\n')
    print('Adding '+str(np.size(added))+' files \n')
    f = 0
    while f < len(added):
        add = ''.join(added[f])   #Convert to string from list
        type = rha2.find_alert_type(path_to_watch, add)
        if type == 'ash':
            out2 = rha2.make_files(path_to_watch, add)
            rha2.run_hysp(out2[0], out2[1])
            rha2.make_disp(out2[0], out2[1], out2[2])
            #rha2.compress_kml(out2[0], out2[1])
            #Check if dispersion figures are generated
            dispfile = rha2.check_for_figs(out2[0],out2[1])
            if dispfile == True:
                original1.append(add)
                print('Dispersion files do exist! \n')
        else:
            original1.append(add)
            print('Not an ash alert: '+add+'\n')
        f += 1

if removed:
    print('Files removed:',', \n'.join(removed),'\n')
    g = 0
    while g < len(removed):
        remov = ''.join(removed[g])
        original1.remove(remov)
        g += 1

if added or removed:
    #Write current list of files that successfully created trajectory figs to original_list2.txt
    print('Updates recorded to file!')
    with open(direct+'/original_list2_disp.txt', 'w') as fis:
                fis.write(json.dumps(original1))
    os.system('mv '+direct+'original_list2_disp.txt '+direct+'original_list_disp.txt')
else:
    print('No updates to '+path_to_watch+' folder!')






