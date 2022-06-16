'''Sample script to easily generate analysis files.
Example command to run script: 

python3 process_files.py kohola 446 448

Calling the script as above processes scan numbers 446-448 for the detector module named "kohola" 
'''

from process_h5 import process_h5
from os import sys

tpc_name = sys.argv[1] #name of TPC module (case sensitive)
first_run = int(sys.argv[2])
last_run = int(sys.argv[3])

for run in range(first_run,last_run+1):
    input_file = 'input_test_files/%s_%s_stop_mode_ext_trigger_scan_interpreted.h5'%(run,tpc_name) #Must include '.h5' file extension. Change input directory as needed
    output_file = '%s_%s_test.root'%(run,tpc_name) #Can use '.root', '.pkl', or '.feather' file extension. Name must include file extension
    print('\n \n Creating ntuple for run %s. Last run is %s \n \n'%(run,last_run))
    
    process_h5(input_file,output_file,save = True) #Uncomment if you haven't explicitly calibrated charge
    #process_h5(input_file,output_file,save = True, calibration_file = 'test_calibration_tables/%s_calibrated.pkl'%(tpc_name)) #Uncomment if you've generated a specific calibration file and add the file name in the 'calibration_file' argument
