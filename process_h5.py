'''
JTS 08/17/2021
Module merges raw hit data and meta_data. Previous scripts discarded many events 
by trying to match event numbers (trigger counts) of metadata tables with hit tables. 
This new script fixes this by keeping all events in hits table; assigning timestamps
to event_numbers that match between hit table and metadata table; and filling all timestamp gaps
in modified hits table by performing linear interpolations between nearest timestamp pairs
'''

import pandas as pd
import numpy as np
import array
import ROOT
from ROOT import TVector3
from tqdm import tqdm
from numba import jit
import os
pd.set_option('mode.chained_assignment', None) #remove pandas copy warning
import warnings
warnings.filterwarnings('ignore') #Used to ignore numba warnings. Comment out when debugging

class process_h5:
    def __init__(self, input_file, output_file, calibration_file = None, save = True, neutron_skim = False):
        self.neutron_skim = neutron_skim
        self.input_file = input_file
        self.output_file = output_file
        #Input detector calibration file
        if calibration_file is not None:
            self.calibration = pd.read_pickle(calibration_file) #Calibration files are stored as pickles
        else:
            self.calibration = pd.read_pickle("test_calibration_tables/test_calibration.pkl")
        
        #Check if filenames are of correct format
        input_extension = os.path.splitext(self.input_file)[1]
        output_extension = os.path.splitext(self.output_file)[1] #grabs file extension from output_file
        if input_extension != '.h5':
            raise NameError("input_file must have a '.h5' extension")
        if output_extension != '.root' and output_extension != '.pkl' and output_extension != '.feather':
            raise NameError("output_file must include '.root', '.pkl', or '.feather' extension")

        #Import hits and metadata tables
        self.raw_hits = pd.read_hdf(input_file, 'Hits') #reads hits into pandas dataframe
        self.meta = pd.read_hdf(input_file, 'meta_data') #reads metadata table

        #Merge hits and metadata tables
        self.hits = self.group_hits_to_event() #grouped hits dataframe
        self.merged_raw = self.merge_hits_and_meta()

        #Process data
        self.processed_data = self.process_merged_data()
        if self.neutron_skim: #Loose selections that keep most neutron events to save storage overhead
            nskim = 'LAPA < 6000 & track_energy > 1 & hitside_row_max == 0 & hitside_col_max == 0'
            +'& hitside_row_min == 0 & hitside_col_min == 0'
            self.processed_data = self.processed_data.query(nskim)
            self.processed_data.index = [i for i in range(0,len(self.processed_data))]

        #Write output data files
        if save:
            if output_extension == '.root': #save as ROOT
                self.write_ROOT_output()
            else: #save as pickle or feather
                self.write_output(output_extension)

    def group_hits_to_event(self): #'unflattens' hits table. Reformats hits table to show entries event by event with hit data stored in arrays and track data stored as numbers
        self.raw_hits['pixel_charge'] = self.raw_hits['tot'].map(self.calibration.set_index('tot_code')['q_per_tot']) #maps tot code to electron charge from calibration file
        self.raw_hits['pixel_energy'] = self.raw_hits['tot'].map(self.calibration.set_index('tot_code')['E_per_tot']) #maps tot code to energy using energy = Q*W/gain ... this function is applied in code generating calibration file
        self.raw_hits['column'] = self.raw_hits['column'] - 1
        self.raw_hits['row'] = self.raw_hits['row'] - 1
        self.raw_hits['x'] = self.raw_hits['column']*250.
        self.raw_hits['y'] = (335 - self.raw_hits['row'])*50. #To make x,y,z a hight handed coordinate system
        try: #If calibration file has a drift speed use it
            self.raw_hits['z'] = self.raw_hits['relative_BCID']*self.calibration['v_drift']
        except:
            self.raw_hits['z'] = self.raw_hits['relative_BCID']*250. #defaults to 250
        print('Grouping hit data into arrays to merge with metadata')
        with tqdm(total=9) as progress:
            col = self.raw_hits.groupby('event_number')['column'].apply(lambda x: np.array(x))
            progress.update(1)
            row = self.raw_hits.groupby('event_number')['row'].apply(lambda x: np.array(x))
            progress.update(1)
            tot = self.raw_hits.groupby('event_number')['tot'].apply(lambda x: np.array(x))
            progress.update(1)
            x = self.raw_hits.groupby('event_number')['x'].apply(lambda x: np.array(x))
            progress.update(1)
            y = self.raw_hits.groupby('event_number')['y'].apply(lambda x: np.array(x))
            progress.update(1)
            z = self.raw_hits.groupby('event_number')['z'].apply(lambda x: np.array(x))
            progress.update(1)
            bc = self.raw_hits.groupby('event_number')['relative_BCID'].apply(lambda x: np.array(x))
            progress.update(1)
            q = self.raw_hits.groupby('event_number')['pixel_charge'].apply(lambda x: np.array(x)) #computes status code for event. OK codes are 8228 or 36, rest should be rejected
            progress.update(1)
            E = self.raw_hits.groupby('event_number')['pixel_energy'].apply(lambda x: np.array(x)) #computes status code for event. OK codes are 8228 or 36, rest should be rejected
            progress.update(1)
            print('DONE!')
        err = self.raw_hits.groupby('event_number')['event_status'].mean() #computes status code for event. OK codes are 8228 or 36, rest should be rejected
        grouped_hits = pd.DataFrame() #populate dataframe with arrays of hit data
        grouped_hits['npoints'] = col.apply(lambda x: len(x)) #number of hits in the event
        grouped_hits['event_number'] = col.index.to_numpy()
        grouped_hits['column'] = col
        grouped_hits['row'] = row
        grouped_hits['BCID'] = bc
        grouped_hits['tot'] = tot
        grouped_hits['x'] = x
        grouped_hits['y'] = y
        grouped_hits['z'] = z
        grouped_hits['evt_status'] = err
        grouped_hits['pixel_charge'] = q
        grouped_hits['pixel_energy'] = E
        grouped_hits.index = [i for i in range(0,len(grouped_hits))] #reindex
        grouped_hits['error_code'] = 0
        index = grouped_hits.query('evt_status != 8228 & evt_status != 36').index.to_numpy() #indices of bad events
        grouped_hits['error_code'][index] = 1 #error_code = 1 means bad event that should be discarded
        return grouped_hits
        
    def merge_hits_and_meta(self): #adds timestamps to grouped hits table
        print('Merging hit data and metadata')
        meta = self.meta.loc[self.meta['event_number'].duplicated() == False] #remove duplicate metadata entries
        merged = self.hits.loc[self.hits['event_number'].duplicated() == False]
        merged['true_ts'] = np.nan #timestamps of hits events matching up with meta events
        ts_array = meta.loc[meta['event_number'].isin(merged['event_number'])]['timestamp_start'].to_numpy()
        true_ts_index = merged.loc[merged['event_number'].isin(meta['event_number'])].index.to_numpy()
        merged['true_ts'][true_ts_index] = ts_array #true ts. May be lots of NANs
        merged['ts'] = merged['true_ts'].interpolate('linear') #linearly interpolated timestamp
        merged = merged.drop(columns = ['true_ts'])
        merged = merged.dropna() #will drop ts entries that weren't interpolated
        merged.index = [i for i in range(0,len(merged))] #reindex after dropping
        cols = [col for col in merged.columns] #list of columns of merged
        cols = cols[-1:] + cols[:-1]
        merged = merged[cols] #move 'ts' to be first column
        merged = merged.query('error_code == 0')
        merged.index = [i for i in range(0,len(merged))]
        return merged

    def process_merged_data(self):

        ### Computations of observables of interest ###

        @jit(nopython=True) #numba for fast computation
        def get_principal_axis(data):
            uu, dd, vv = np.linalg.svd(data-np.array([data[:,0].mean(),data[:,1].mean(),data[:,2].mean()]))
            projection = (data @ vv.T).T[0]
            return projection.max() - projection.min(), vv[0] #returns track length and principal axis vector

        ##Transform carefully to head and tail direction

        def get_angles_and_charge_fractions(data,vec,q):
            Tvec = TVector3(vec[0],vec[1],vec[2])
            theta = Tvec.Theta() #RADIANS
            phi = Tvec.Phi() #RADIANS
            if np.cos(theta) < 0: #restrict vector so head points up always
                vec = -1 * Tvec
                theta = Tvec.Theta()
                phi = Tvec.Phi()
            v = np.array([Tvec.x(), Tvec.y(), Tvec.z()]) #convert back to numpy
            projection = data @ v.T
            midp = 0.5*float(projection.max()+projection.min())
            uc = 0 #upper half charge
            lc = 0 #lower half charge
            for i,val in enumerate(projection):
                if val > midp:
                    uc += q[i]
                elif val < midp:
                    lc += q[i]
                elif val == midp and i%2 == 0:
                    uc += q[i]
                elif val == midp and i%2 != 0:
                    lc += q[i]
            upper_charge_fraction = uc/(uc+lc)
                
            return theta, phi, upper_charge_fraction, uc, lc
    
        @jit(nopython=True)
        def SDCD(data): #standard deviation of charge distribution
            return np.sqrt(np.diag(data @ data.T).mean())
    
        @jit(nopython=True)
        def wSDCD(wdata,q): #charge weighted standard deviation of charge distribution
            return np.sqrt(np.diag(wdata @ wdata.T).sum()/(q.sum()))

        @jit(nopython=True)
        def compute_y_track_hat(zhat,xhat): #"track coordinate" y unit vector
            return np.cross(zhat,xhat)

        @jit(nopython=True)
        def compute_z_track_hat(xhat,yhat): #"track coordinate" z unit vector
            return np.cross(xhat,yhat)

        @jit(nopython=True)
        def compute_track(data,tuv): #computes track coordinate system; tuv is track unit vector
            return data @ tuv

        def ChargeUnif(data): #std dev of distribution of mean distances between each charge and all other charges
            a = np.linalg.norm(data - data[:,None], axis=-1)
            return np.std([a[i].mean() for i in range(0,len(data))])

        # initial dataframe computations
        
        data = self.merged_raw
        data['sum_tot'] = data['tot'].apply(lambda x: x.sum()) #total TOT in event
        data['sat_frac'] = data['tot'].apply(lambda x: len(np.where(x == 13)[0])/len(x)) #computes fraction of pixels with TOT = 13
        data['track_charge'] = data['pixel_charge'].apply(lambda x: x.sum()) #sums all hit charges to get total charge of event
        data['track_energy'] = data['pixel_energy'].apply(lambda x: x.sum()) #sums all hit energies to get total energy of event
        if self.neutron_skim: #remove lowest energy events to speed up processing
            data = data.query('track_energy > 1')
            data.index = [i for i in range(0,len(data))]
        data['hitside_row_min'] = data['row'].apply(lambda x: 0 if len(np.where(x == 0)[0]) == 0 else 1) #1 if bottom row edge pixel is hit
        data['hitside_row_max'] = data['row'].apply(lambda x: 0 if len(np.where(x == 335)[0]) == 0 else 1) #1 if top row edge pixel is hit
        data['hitside_col_min'] = data['column'].apply(lambda x: 0 if len(np.where(x == 0)[0]) == 0 else 1) #1 if left column edge pixel is hit
        data['hitside_col_max'] = data['column'].apply(lambda x: 0 if len(np.where(x == 79)[0]) == 0 else 1) #1 if right column edge pixel is hit
        data['x_center'] = data['x'] - data['x'].apply(lambda x: x.mean()) #center tracks for track coordinate computations
        data['y_center'] = data['y'] - data['y'].apply(lambda x: x.mean())
        data['z_center'] = data['z'] - data['z'].apply(lambda x: x.mean())

        # lists to fill

        ls = [] #lengths
        xt_hats = [] #principal axis vectors, direction of x_track_hat
        yt_hats = []
        zt_hats = []
        xtracks = [] #coordinates in track coordinate system
        ytracks = []
        ztracks = []
        stds = [] #std deviations of charge distribution
        wstds = [] #charge weighted std deviations of charge distribution
        thetas = [] #zenith angles
        phis = [] #azimuthal angles with respect to readout plane
        ucfs = [] #Fractions of total track charge on upper half of track
        ucs = []
        lcs = []
        chargeunif = []
        
        print("Computing e rejection discriminants...")
        zhat = np.array([0,0,1],dtype='float32') #zhat in detector coordinates
        
        for i in tqdm(range(0,len(data))):
            track = np.concatenate([[data['x'][i].T,data['y'][i].T,data['z'][i].T]]).T
            ctrack = np.concatenate([[data['x_center'][i].T,data['y_center'][i].T,data['z_center'][i].T]]).T #centered track
            wtrack = np.concatenate([[data['pixel_charge'][i]*data['x_center'][i].T,data['pixel_charge'][i]*data['y_center'][i].T,data['pixel_charge'][i]*data['z_center'][i].T]]).T # for wSDCD
            q = data['pixel_charge'][i]
            l, xt_hat = get_principal_axis(track) #length and principal axis vector
            ls.append(l)
            xt_hats.append(xt_hat) #principal axis direction
            yt_hat = compute_y_track_hat(zhat,xt_hat)
            ynorm = np.linalg.norm(yt_hat)
            yt_hat = yt_hat/ynorm #normalize
            yt_hats.append(yt_hat)
            zt_hat = compute_z_track_hat(xt_hat,yt_hat)
            zt_hats.append(zt_hat)
            xtrack = compute_track(ctrack,xt_hat)
            ytrack = compute_track(ctrack,yt_hat)
            ztrack = compute_track(ctrack,zt_hat)
            xtracks.append(xtrack)
            ytracks.append(ytrack)
            ztracks.append(ztrack)
            stds.append(SDCD(ctrack))
            wstds.append(wSDCD(wtrack,q))
            chargeunif.append(ChargeUnif(track))
            th, ph, ucf, uc, lc = get_angles_and_charge_fractions(track,xt_hat,q) #theta, phi, upper charge fraction, upper charge, lower charge, respectively
            thetas.append(th)
            phis.append(ph)
            ucfs.append(ucf)
            ucs.append(uc)
            lcs.append(lc)
        data['x_track'] = xtracks
        data['y_track'] = ytracks
        data['z_track'] = ztracks
        data['theta'] = thetas
        data['phi'] = phis
        data['upper_charge_fraction'] = ucfs
        data['upper_charge'] = ucs
        data['lower_charge'] = lcs
        data['LAPA'] = ls
        data['SDCD'] = stds
        data['wSDCD'] = wstds
        data['CylThick'] = (data['y_track']**2+data['z_track']**2).apply(lambda x: x.sum())
        data['ChargeUnif'] = chargeunif
        #data['PrincipalAxis'] = xt_hats #direction along PA
        print('DONE!')        
        return data

    # Routine for making ROOT ntuple

    def write_ROOT_output(self):
        data = self.processed_data
        data['raw_event_number'] = data['event_number']
        data['event_number'] = data.index
        keys = [val for val in data.columns]
        keys = keys[-1:] + keys[:-1] #reorder columns to move raw_event_number in front
        data = data[keys] #reorder columns to move raw_event_number in front
        output = ROOT.TFile(self.output_file, 'recreate')
        tout = ROOT.TTree('data','data')
        branches = {}
        root_data={}        
        for key in keys:
            if data[key].dtype == "O": #Determines the size of an array in a dataframe to be pushed to ntuple
                root_data[key]=array.array('d',[0 for j in range(0,100000)])
                branches[key]=tout.Branch("%s"%(key), root_data[key], "%s[npoints]/D"%(key))
            elif key == 'npoints':
                root_data[key]=array.array('i',[0])
                branches[key]=tout.Branch("%s"%(key), root_data[key], "%s/I"%(key))
            else:
                root_data[key]=array.array('d',[0])
                branches[key]=tout.Branch("%s"%(key), root_data[key], "%s/D"%(key))
        print('Filling ROOT Tree')
        for j in tqdm(range(0,len(data))):
            root_data['npoints'][0] = data['npoints'].to_numpy()[j].astype(int)
            for key in keys:
                if data[key].dtype == "O":
                    for i in range(0,root_data['npoints'][0]):
                        root_data[key][i]=data[key][j][i]
                elif key != 'npoints':
                    root_data[key][0]=data[key][j]
            tout.Fill()
        print('DONE!')
        output.Write()
        output.Close()

    def write_output(self,output_extension):
        data = self.processed_data
        data['raw_event_number'] = data['event_number']
        data['event_number'] = data.index
        keys = [val for val in data.columns]
        keys = keys[-1:] + keys[:-1] #reorder columns to move raw_event_number in front
        data = data[keys] #reorder columns to move raw_event_number in front
        if output_extension == '.pkl':
            data.to_pickle(self.output_file)
        else:
            data.to_feather(self.output_file)

    if __name__ == '__main__': #Recommend using this as a module rather than standalone code
        process_h5()
