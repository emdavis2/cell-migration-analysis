from functions.PRW_model_functions import *
from functions.PRWpolaritybias_model_functions import *
from functions.acf_functions import *

import numpy as np 

def run_sim_get_err(poslagaverage_data, run_sim_fn, min_track_length):
    #np.random.seed(20)

    data_sim = run_sim_fn

    poslagaverage = np.zeros(300)
    all_ac = []
    for df in data_sim:
        track=df
        combined = make_comb_df(track['vx'].to_list()[2:min_track_length-2],track['vy'].to_list()[2:min_track_length-2])
        combined = combined.dropna()
        poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined)

        #remove nans here
        poslagsmean[np.isnan(poslagsmean)] = 0
        all_ac.append(poslagsmean)
        poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
    poslagaverage /= len(data_sim) #Nposlagtotal 

    poslagaverage_sim = poslagaverage[0:min_track_length-4]

    acf_vel_err = np.sum(np.abs(poslagaverage_data - poslagaverage_sim))

    return acf_vel_err


def perform_gridsearch_2params(poslagaverage_data, std_dev_w_vals, std_dev_theta_vals, Nwalkers, dt, time, min_track_length):
    min_err = 100
    index_w = 0
    index_theta = 0

    for ind_w, std_dev_w in enumerate(std_dev_w_vals):
        for ind_t,std_dev_theta in enumerate(std_dev_theta_vals):
            err = run_sim_get_err(poslagaverage_data, run_PRWpolaritybias_sim(Nwalkers, dt, time, std_dev_w, std_dev_theta), min_track_length)
            if err < min_err:
                min_err = err
                index_w = ind_w
                index_theta = ind_t
            else:
                continue
    
    return min_err, std_dev_w_vals[index_w], std_dev_theta_vals[index_theta]


def perform_gridsearch_1param(poslagaverage_data, std_dev_theta_vals, Nwalkers, dt, time, min_track_length):
    min_err = 100
    index_theta = 0

    for ind_t, std_dev_theta in enumerate(std_dev_theta_vals):
        err = run_sim_get_err(poslagaverage_data, run_PRW_sim(Nwalkers, dt, time, std_dev_theta), min_track_length)
        if err < min_err:
            min_err = err
            index_theta = ind_t
        else:
            continue
    
    return min_err, std_dev_theta_vals[index_theta]