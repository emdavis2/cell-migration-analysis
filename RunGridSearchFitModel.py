from functions.acf_functions import *
from functions.compile_data_tracks_function import *
from functions.model_fitting_functions import *
from functions.PRW_model_functions import *
from functions.PRWpolaritybias_model_functions import *
from functions.langevin_PRW_functions import *

import sys
import numpy as np
import pandas as pd

treatment = str(sys.argv[1])

min_track_length = int(sys.argv[2])

region = str(sys.argv[3])

time = int(sys.argv[4]) #5

dt = float(sys.argv[5]) #0.1667

Nwalkers = int(sys.argv[6]) #113

model_type = str(sys.argv[7]) #PRW or PRW_PB or LPRW

err_fn = str(sys.argv[8]) #vel_acf or MSD

param1_name = str(sys.argv[9])

param1_start = float(sys.argv[10])

param1_stop = float(sys.argv[11])

param1_num = int(sys.argv[12])

param2_name = str(sys.argv[13])

param2_start = float(sys.argv[14])

param2_stop = float(sys.argv[15])

param2_num = int(sys.argv[16])

tracks_region, tracks_geo_region, region_cells, region_endpointcells = compile_data_tracks(treatment, min_track_length, region)


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

#MSD for data
MSD_data = calc_MSD(tracks_region, min_track_length)

if err_fn == 'vel_acf':
    data_for_fit = poslagaverage_data
elif err_fn == 'MSD':
    data_for_fit = MSD_data

#perform grid search
if model_type == 'PRW':
    std_dev_theta_vals = np.linspace(0.4, 1.5, 20)
    min_err, std_dev_theta = perform_gridsearch_1param(data_for_fit, model_type, err_fn, std_dev_theta_vals, Nwalkers, dt, time, min_track_length)
    with open(r'model/model_params_{}_{}_{}.txt'.format(region, model_type, err_fn), 'w') as f:
        f.write('min_err={}'.format(str(min_err)))
        f.write('\n')
        f.write('std_dev_theta={}'.format(str(std_dev_theta)))
elif model_type == 'PRW_PB' or 'LPRW':
    #std_dev_w_vals = np.linspace(0.2, 0.9, 20)
    #std_dev_theta_vals = np.linspace(0.9, 1.5, 20)
    param1_vals = np.linspace(param1_start,param1_stop,param1_num)
    param2_vals = np.linspace(param2_start,param2_stop,param2_num)
    min_err, param1, param2 = perform_gridsearch_2params(data_for_fit, model_type, err_fn, param1_vals, param2_vals, Nwalkers, dt, time, min_track_length)
    with open(r'model/model_params_{}_{}_{}.txt'.format(region, model_type, err_fn), 'w') as f:
        f.write('min_err={}'.format(str(min_err)))
        f.write('\n')
        f.write('{}={}'.format(str(param1_name), str(param1)))
        f.write('\n')
        f.write('{}={}'.format(str(param2_name), str(param2)))

