from functions.PRW_model_functions import *
from functions.PRWpolaritybias_model_functions import *
from functions.weighted_PRW_model_functions import *
from functions.langevin_PRW_functions import *
from functions.acf_functions import *
from functions.msd_functions import *

import numpy as np 
import pandas as pd

def make_comb_df(vx, vy):
    d = {'vx': vx, 'vy': vy}
    vel_df=pd.DataFrame(data=d)
    combined = pd.concat([vel_df[['vx','vy']].reset_index(drop=True), vel_df[['vx','vy']].reset_index(drop=True)], axis = 1 )
    return combined

def run_sim_get_velacf_err(poslagaverage_data, run_sim_fn, min_track_length):
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

def run_sim_get_MSD_err(MSD_data, run_sim_fn, min_track_length):
  data_sim = run_sim_fn

  MSD_sim = calc_MSD_sim(data_sim, min_track_length)

  MSD_err = np.sum(np.abs(MSD_data - MSD_sim))

  return MSD_err

def perform_gridsearch_3params(data_for_fit, run_sim_fn, run_sim_err_fn, param1_vals, param2_vals, param3_vals, Nwalkers, dt, time, min_track_length):
    dispatcher_err = {'vel_acf': run_sim_get_velacf_err, 'MSD': run_sim_get_MSD_err}
    dispatcher_sim = {'weighted_PRW': run_weighted_PRW_sim}
    min_err = 1000000
    index_p1 = 0
    index_p2 = 0
    index_p3 = 0

    for ind_p1, param1 in enumerate(param1_vals):
        for ind_p2, param2 in enumerate(param2_vals):
            for ind_p3, param3 in enumerate(param3_vals):
                err = dispatcher_err[run_sim_err_fn](data_for_fit, dispatcher_sim[run_sim_fn](Nwalkers, dt, time, param1, param2, param3), min_track_length)
                if err < min_err:
                    min_err = err
                    index_p1 = ind_p1
                    index_p2 = ind_p2
                    index_p3 = ind_p3
                else:
                    continue
    
    return min_err, param1_vals[index_p1], param2_vals[index_p2], param3_vals[index_p3]

def perform_gridsearch_2params(data_for_fit, run_sim_fn, run_sim_err_fn, param1_vals, param2_vals, Nwalkers, dt, time, min_track_length):
    dispatcher_err = {'vel_acf': run_sim_get_velacf_err, 'MSD': run_sim_get_MSD_err}
    dispatcher_sim = {'PRW_PB': run_PRWpolaritybias_sim, 'LPRW': run_PRW_langevin_sim}
    min_err = 1000000
    index_p1 = 0
    index_p2 = 0

    for ind_p1, param1 in enumerate(param1_vals):
        for ind_p2, param2 in enumerate(param2_vals):
            err = dispatcher_err[run_sim_err_fn](data_for_fit, dispatcher_sim[run_sim_fn](Nwalkers, dt, time, param1, param2), min_track_length)
            if err < min_err:
                min_err = err
                index_p1 = ind_p1
                index_p2 = ind_p2
            else:
                continue
    
    return min_err, param1_vals[index_p1], param2_vals[index_p2]


def perform_gridsearch_1param(data_for_fit, run_sim_fn, run_sim_err_fn, std_dev_theta_vals, Nwalkers, dt, time, min_track_length):
    dispatcher_err = {'vel_acf': run_sim_get_velacf_err, 'MSD': run_sim_get_MSD_err}
    dispatcher_sim = {'PRW':run_PRW_sim}
    min_err = 1000000
    index_theta = 0

    for ind_t, std_dev_theta in enumerate(std_dev_theta_vals):
        err = dispatcher_err[run_sim_err_fn](data_for_fit, dispatcher_sim[run_sim_fn](Nwalkers, dt, time, std_dev_theta), min_track_length)
        if err < min_err:
            min_err = err
            index_theta = ind_t
        else:
            continue
    
    return min_err, std_dev_theta_vals[index_theta]