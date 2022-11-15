from functions.PRW_model_functions import *
from functions.PRWpolaritybias_model_functions import *
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

def run_sim_get_velacf_err(tracks_region, tracks_geo_region, run_sim_fn, min_track_length, walker):
    #autocorrelation velocity for data
    track_data = tracks_geo_region[walker]
    combined_data = pd.concat([track_data[['vx','vy']].reset_index(drop=True), track_data[['vx','vy']].reset_index(drop=True)], axis = 1 )
    poslagsmean_data, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined_data, min_track_length)
    poslagsmean_data[np.isnan(poslagsmean_data)] = 0

    track = run_sim_fn[0]

    combined = make_comb_df(track['vx'].to_list()[2:min_track_length-2],track['vy'].to_list()[2:min_track_length-2])
    combined = combined.dropna()
    poslagsmean_sim, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

    #remove nans here
    poslagsmean_sim[np.isnan(poslagsmean_sim)] = 0

    acf_vel_err = np.sum(np.abs(poslagsmean_data - poslagsmean_sim))

    return acf_vel_err

def run_sim_get_MSD_err(tracks_region, tracks_geo_region, run_sim_fn, min_track_length, walker):

  data_pos = {'x': tracks_region[walker][:,0],'y':tracks_region[walker][:,1]}
  datapos_df = pd.DataFrame(data=data_pos)
  MSD_data_onecell = calc_MSD_onecell(datapos_df, min_track_length)

  data_sim = run_sim_fn[0]

  MSD_sim = calc_MSD_onecell(data_sim, min_track_length)

  MSD_err = np.sum(np.abs(MSD_data_onecell - MSD_sim))

  return MSD_err


def perform_gridsearch_2params(tracks_region, tracks_geo_region, run_sim_fn, run_sim_err_fn, param1_vals, param2_vals, Nwalkers, dt, time, min_track_length):
    dispatcher_err = {'vel_acf': run_sim_get_velacf_err, 'MSD': run_sim_get_MSD_err}
    dispatcher_sim = {'PRW_PB': run_PRWpolaritybias_sim, 'LPRW': run_PRW_langevin_sim}
    tot_err = 0
    index_p1_allcells = []
    index_p2_allcells = []
    for walker in range(Nwalkers):
        min_err = 100
        index_p1 = 0
        index_p2 = 0

        for ind_p1, param1 in enumerate(param1_vals):
            for ind_p2,param2 in enumerate(param2_vals):
                err = dispatcher_err[run_sim_err_fn](tracks_region, tracks_geo_region, dispatcher_sim[run_sim_fn](1, dt, time, param1, param2), min_track_length, walker)
                if err < min_err:
                    min_err = err
                    index_p1 = ind_p1
                    index_p2 = ind_p2
                else:
                    continue
        index_p1_allcells.append(index_p1)
        index_p2_allcells.append(index_p2)
        tot_err += min_err

    param1_list = [param1_vals[i] for i in index_p1_allcells]
    param2_list = [param2_vals[i] for i in index_p2_allcells]
    
    return tot_err, param1_list, param2_list


def perform_gridsearch_1param(tracks_region, tracks_geo_region, run_sim_fn, run_sim_err_fn, std_dev_theta_vals, Nwalkers, dt, time, min_track_length):
    dispatcher_err = {'vel_acf': run_sim_get_velacf_err, 'MSD': run_sim_get_MSD_err}
    dispatcher_sim = {'PRW':run_PRW_sim}
    tot_err = 0
    index_theta_allcells = []
    for walker in range(Nwalkers):
        min_err = 100
        index_theta = 0

        for ind_t, std_dev_theta in enumerate(std_dev_theta_vals):
            err = dispatcher_err[run_sim_err_fn](tracks_region, tracks_geo_region, dispatcher_sim[run_sim_fn](1, dt, time, std_dev_theta), min_track_length, walker)
            if err < min_err:
                min_err = err
                index_theta = ind_t
            else:
                continue
        index_theta_allcells.append(index_theta)
        tot_err += min_err

    theta_list = [std_dev_theta_vals[i] for i in index_theta_allcells]
    
    return tot_err, theta_list