
**Note: [ROOT](https://root.cern/install/) is required to run these data processing scripts. I would recommend creating a new environment in anaconda and installing ROOT in that environment.**

# tpc_data_processing
This repository processes FE-I4 pixel chip data recorded using the pyBAR software package. The repository was built for the BEAST TPC detectors and includes the following:  

1. process_h5.py: A module that reads pyBAR interpreted hdf5 files and processes them into either ROOT, pickle, or Apache Arrow feather files. The processed files include event level information such as: raw pixel-hit data, and timestamps. This module also uses a singular value decomposition to identify the principal axis of the track and computes several shape variables that are of interest to various analysis applications  
2. process_files.py: A simple script that uses the process_h5 module to generate output files given numbered input files
3. calibrate_charge.py: A script that takes the output of pyBAR's calibrate_tot.py script and generates a lookup table to map pixel level charge units to physical charge units.If detector gain has been calibrated, this script can also calibrate the energy scale of the TPCs  
 
 
In general the user would generate all relevant files to interface with these scripts using pyBAR, however sample .h5 data files and calibration files are included to test these scripts and generate output files for analysis.
