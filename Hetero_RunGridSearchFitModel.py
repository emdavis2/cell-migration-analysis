from functions.acf_functions import *
from functions.compile_data_tracks_function import *
from functions.hetero_model_fitting_functions import *

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

#perform grid search
if model_type == 'PRW':
    std_dev_theta_vals = np.linspace(0.1, 2, 10)
    tot_err, std_dev_theta = perform_gridsearch_1param(tracks_region, tracks_geo_region, model_type, err_fn, std_dev_theta_vals, Nwalkers, dt, time, min_track_length)
    with open(r'hetero_model/model_params_{}_{}_{}.txt'.format(region, model_type, err_fn), 'w') as f:
        f.write('min_err={}'.format(str(tot_err)))
        f.write('\n')
        f.write('std_dev_theta={}'.format(str(std_dev_theta)))
elif model_type == 'PRW_PB' or 'LPRW':
    #std_dev_w_vals = np.linspace(0.1, 2, 10)
    #std_dev_theta_vals = np.linspace(0.1, 2, 10)
    param1_vals = np.linspace(param1_start,param1_stop,param1_num)
    param2_vals = np.linspace(param2_start,param2_stop,param2_num)
    tot_err, param1, param2 = perform_gridsearch_2params(tracks_region, tracks_geo_region, model_type, err_fn, param1_vals, param2_vals, Nwalkers, dt, time, min_track_length)
    with open(r'hetero_model/model_params_{}_{}_{}.txt'.format(region, model_type, err_fn), 'w') as f:
        f.write('min_err={}'.format(str(tot_err)))
        f.write('\n')
        f.write('{}={}'.format(str(param1_name), str(param1)))
        f.write('\n')
        f.write('{}={}'.format(str(param2_name), str(param2)))

