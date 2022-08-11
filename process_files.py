'''Sample script to easily generate analysis files.
Example command to run script: 

python3 process_files.py kohola 446 448 ROOT

Calling the script as above processes scan numbers 446-448 for the detector module named "kohola" 
'''

from process_h5 import process_h5
import os
from os import sys

tpc_name = sys.argv[1] #name of TPC module (case sensitive)
first_run = int(sys.argv[2])
last_run = int(sys.argv[3])
output_extension = sys.argv[4] # Use either "root", "feather", or "pickle"

if len(sys.argv) != 5 or output_extension.lower() not in ['root', 'pickle', 'feather']:
    raise ValueError("\n\nScript must be called as python3 process_files.py <tpc_name> <first_run> <last_run> <output_file_extension_type>. Acceptable arguments for <output_file_extension_type> (case not sensitive): 'ROOT', 'pickle', or 'feather'.\n\nAs an example:   python3 process_files.py kohola 446 448 feather\n")

for run in range(first_run,last_run+1):
    
    input_file = 'input_test_files/%s_%s_stop_mode_ext_trigger_scan_interpreted.h5'%(run,tpc_name) #Must include '.h5' file extension. Change input directory as needed

    if os.path.splitext(input_file)[1] != '.h5':
        raise ValueError("Input file must be .h5 . By default we read in *_stop_mode_ext_trigger_scan_interpreted.h5, so look for files with this naming convention!")

    if output_extension.lower() == 'root':
        output_file = 'output_analysis_files/ROOT/%s_%s_test.root'%(run,tpc_name)
    elif output_extension.lower() == 'pickle':
        output_file = 'output_analysis_files/pickle/%s_%s_test.pkl'%(run,tpc_name)
    elif output_extension.lower() == 'feather':
        output_file = 'output_analysis_files/feather/%s_%s_test.feather'%(run,tpc_name)
    
    print('\n \n Creating ntuple for run %s. Last run is %s \n \n'%(run,last_run))
    
    process_h5(input_file,output_file,save = True) #Uncomment if you haven't explicitly calibrated charge

    ### Comment previous line and uncomment the following line if you've generated a specific calibration file and add the file name in the 'calibration_file' argument
    #process_h5(input_file,output_file,save = True, calibration_file = 'test_calibration_tables/%s_calibrated.pkl'%(tpc_name)) 
