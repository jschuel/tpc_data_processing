''' Reads a tot calibration h5 file and extracts the mean tot code for a set of plsrdac charges as given in the file
 As an example, this script could be run with python3 calibrate_tpc.py kohola 15
'''
import sys
import h5py
import numpy as np
import tables
import ROOT
import pandas as pd
import array

tpc = sys.argv[1]
run_num = sys.argv[2]

def create_lookup_table(tpc, scan_number):
    filename = 'test_calibration_input/%s_%s_tot_calibration_interpreted.h5' % (scan_number,tpc)
    f = h5py.File(filename, 'r')
    df = pd.read_hdf(filename, 'meta_data') #metadata table can be directly read into pandas dataframe
    mean_tot = f['HistMeanTot'][:] #(336 x 80 x #plsrdacbins) ndarray
    mean_tot[np.isnan(mean_tot)]=0 #masks nan with 0
    means = np.mean(np.mean(mean_tot,axis=0),axis =0) #inner mean is over rows, outer mean over columns
    err1 = np.sqrt(np.sum(np.std(mean_tot,axis=0)**2,axis=0))/np.sqrt(100) #std_error over # injections
    err2 = np.std(np.mean(mean_tot,axis=0),axis =0) #standard deviation over row mean
    errs = np.sqrt(err1**2+err2**2)
    lt = pd.DataFrame()
    lt['plsrdac'] = df['PlsrDAC'].unique()
    lt['tot_mean'] = means
    lt['tot_err'] = errs
    q_per_tot = charge_per_tot(lt['tot_mean'], lt['tot_err'], lt['plsrdac'], tpc)

    return lt, q_per_tot

def charge_per_tot(mean_tot, err_tot, plsrdac_array, tpc):
    qe = 52 #conversion factor from digital units of charge to physical charge
    W = 34.45 #He:CO2 work function
    electron_charge = qe * plsrdac_array
    tot_code = array.array('i', [i for i in range(0,14)])
    #Create ROOT TGraph container for bilinear interpolation
    gr = ROOT.TGraphErrors(len(mean_tot), array.array('d', mean_tot), array.array('d', electron_charge), array.array('d', err_tot), array.array('d', [0 for i in range(0,len(mean_tot))]))
    df = pd.DataFrame()
    df['tot_code'] = tot_code
    df['q_per_tot'] = [gr.Eval(tot_code[i]) for i in range(0,14)] #Interpolation to map integer time over threshold (TOT) codes to physical charge
    df['E_per_tot'] = df['q_per_tot']*W/909*1e-3 #Convert to energy given a calibrated gain of 909
    return df

table, charge = create_lookup_table("%s"%(tpc), run_num)
charge.to_pickle("test_calibration_tables/%s_calibrated.pkl"%(tpc)) #save charge calibration table as pickle
