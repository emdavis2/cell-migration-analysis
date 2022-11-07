from functions.PRW_model_functions import *
from functions.PRWpolaritybias_model_functions import *
from functions.acf_functions import *
from functions.msd_functions import *

import numpy as np 
import pandas as pd

def make_comb_df(vx, vy):
    d = {'vx': vx, 'vy': vy}
    vel_df=pd.DataFrame(data=d)
    combined = pd.concat([vel_df[['vx','vy']].reset_index(drop=True), vel_df[['vx','vy']].reset_index(drop=True)], axis = 1 )
    return combined

def run_sim_get_velacf_err(tracks_region, tracks_geo_region, run_sim_fn, min_track_length):

    #autocorrelation velocity for data
    poslagaverage_data = np.zeros(300)
    for df in tracks_geo_region:
        track=df
        combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
        poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)
        #remove nans here
        poslagsmean[np.isnan(poslagsmean)] = 0
        poslagaverage_data[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
    poslagaverage_data /= len(tracks_geo_region) #Nposlagtotal 
    poslagaverage_data = poslagaverage_data[0:min_track_length-4]

    data_sim = run_sim_fn

    poslagaverage = np.zeros(300)
    for df in data_sim:
        track=df
        combined = make_comb_df(track['vx'].to_list()[2:min_track_length-2],track['vy'].to_list()[2:min_track_length-2])
        combined = combined.dropna()
        poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

        #remove nans here
        poslagsmean[np.isnan(poslagsmean)] = 0
        poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
    poslagaverage /= len(data_sim) #Nposlagtotal 
    poslagaverage_sim = poslagaverage[0:min_track_length-4]

    acf_vel_err = np.sum(np.abs(poslagaverage_data - poslagaverage_sim))

    return acf_vel_err

def run_sim_get_MSD_err(tracks_region, tracks_geo_region, run_sim_fn, min_track_length):

  MSD_data = calc_MSD(tracks_region, min_track_length)

  data_sim = run_sim_fn

  MSD_sim = calc_MSD_sim(data_sim, min_track_length)

  MSD_err = np.sum(np.abs(MSD_data - MSD_sim))

  return MSD_err


def perform_gridsearch_2params(tracks_region, tracks_geo_region, run_sim_err_fn, std_dev_w_vals, std_dev_theta_vals, Nwalkers, dt, time, min_track_length):
    dispatcher = {'vel_acf': run_sim_get_velacf_err, 'MSD': run_sim_get_MSD_err}
    min_err = 100
    index_w = 0
    index_theta = 0

    for ind_w, std_dev_w in enumerate(std_dev_w_vals):
        for ind_t,std_dev_theta in enumerate(std_dev_theta_vals):
            err = dispatcher[run_sim_err_fn](tracks_region, tracks_geo_region, run_PRWpolaritybias_sim(Nwalkers, dt, time, std_dev_w, std_dev_theta), min_track_length)
            if err < min_err:
                min_err = err
                index_w = ind_w
                index_theta = ind_t
            else:
                continue
    
    return min_err, std_dev_w_vals[index_w], std_dev_theta_vals[index_theta]


def perform_gridsearch_1param(tracks_region, tracks_geo_region, run_sim_err_fn, std_dev_theta_vals, Nwalkers, dt, time, min_track_length):
    dispatcher = {'vel_acf': run_sim_get_velacf_err, 'MSD': run_sim_get_MSD_err}
    min_err = 100
    index_theta = 0

    for ind_t, std_dev_theta in enumerate(std_dev_theta_vals):
        err = dispatcher[run_sim_err_fn](tracks_region, tracks_geo_region, run_PRW_sim(Nwalkers, dt, time, std_dev_theta), min_track_length)
        if err < min_err:
            min_err = err
            index_theta = ind_t
        else:
            continue
    
    return min_err, std_dev_theta_vals[index_theta]