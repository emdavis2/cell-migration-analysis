from functions.acf_functions import *
from functions.compile_data_tracks_function import *
from functions.PRW_model_functions import *
from functions.PRWpolaritybias_model_functions import *
from functions.model_fitting_functions import *

import sys
import numpy as np
import matplotlib.pyplot as plt

treatment = str(sys.argv[1])

min_track_length = int(sys.argv[2])

region = str(sys.argv[3])

time = int(sys.argv[4]) #5

dt = float(sys.argv[5]) #0.1667

Nwalkers = int(sys.argv[6]) #113

PRW_params = str(sys.argv[7])

PRWPB_params = str(sys.argv[8])

#get params from model fitting
PRW_params_open = open(PRW_params, 'r')
PRW_params_readlines = PRW_params_open.readlines()
PRW_err = float(PRW_params_readlines[0].split('=')[1])
PRW_theta_std_dev = float(PRW_params_readlines[1].split('=')[1])

PRWPB_params_open = open(PRWPB_params, 'r')
PRWPB_params_readlines = PRWPB_params_open.readlines()
PRWPB_err = float(PRWPB_params_readlines[0].split('=')[1])
PRWPB_w_std_dev = float(PRWPB_params_readlines[1].split('=')[1])
PRWPB_theta_std_dev = float(PRWPB_params_readlines[2].split('=')[1])

tracks_region, tracks_geo_region, region_cells, region_endpointcells = compile_data_tracks(treatment, min_track_length, region)

#autocorrelation velocity for data
poslagaverage_data = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage_data[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
poslagaverage_data /= len(tracks_geo_region) #Nposlagtotal 
poslagaverage_data = poslagaverage_data[0:min_track_length-4]

#autocorrelation velocity for PRW model 
data_PRWsim = run_PRW_sim(Nwalkers, dt, time, PRW_theta_std_dev)

poslagaverage_PRWsim = np.zeros(300)
all_ac = []
for df in data_PRWsim:
    track=df
    combined = make_comb_df(track['vx'].to_list()[2:min_track_length-2],track['vy'].to_list()[2:min_track_length-2])
    combined = combined.dropna()
    poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

    #remove nans here
    poslagsmean[np.isnan(poslagsmean)] = 0
    all_ac.append(poslagsmean)
    poslagaverage_PRWsim[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
poslagaverage_PRWsim /= len(data_PRWsim) #Nposlagtotal 

poslagaverage_PRWsim = poslagaverage_PRWsim[0:min_track_length-4]

#autocorrelation velocity for PRW polarity bias model
data_PRWPBsim = run_PRWpolaritybias_sim(Nwalkers, dt, time, PRWPB_w_std_dev, PRWPB_theta_std_dev)

poslagaverage_PRWPBsim = np.zeros(300)
all_ac = []
for df in data_PRWPBsim:
    track=df
    combined = make_comb_df(track['vx'].to_list()[2:min_track_length-2],track['vy'].to_list()[2:min_track_length-2])
    combined = combined.dropna()
    poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

    #remove nans here
    poslagsmean[np.isnan(poslagsmean)] = 0
    all_ac.append(poslagsmean)
    poslagaverage_PRWPBsim[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
poslagaverage_PRWPBsim /= len(data_PRWPBsim) #Nposlagtotal 

poslagaverage_PRWPBsim = poslagaverage_PRWPBsim[0:min_track_length-4]

#Compare acf velocity plots
plt.plot(poslagaverage_data,label = "Data")
plt.plot(poslagaverage_PRWPBsim,label = "PRW Polarity Bias, error={}".format(round(PRWPB_err,3)))
plt.plot(poslagaverage_PRWsim,label = "PRW, error={}".format(round(PRW_err,3)))
plt.hlines(y=0,xmin=0,xmax=min_track_length)
plt.xlim(0,min_track_length-2)
plt.ylim(-1,1)
plt.xlabel("time lag")
plt.title(" Autocorrelation velocity {}".format(region))
plt.legend()
plt.savefig('figures/model/acf_velocity_modelcomaparion_{}'.format(region))
plt.clf()

#compare dx and dy with cdf
dx_PRWsim = []
dy_PRWsim = []
for i in range(len(data_PRWsim)):
  dx_PRWsim.append(np.diff(np.array(data_PRWsim[i]['x'].dropna().tolist())))
  dy_PRWsim.append(np.diff(np.array(data_PRWsim[i]['y'].dropna().tolist())))

dx_PRWsim = np.concatenate(dx_PRWsim).ravel()
dy_PRWsim = np.concatenate(dy_PRWsim).ravel()

dx_dy_PRWsim = np.concatenate((dx_PRWsim,dy_PRWsim))

dx_PRWPBsim = []
dy_PRWPBsim = []
for i in range(len(data_PRWPBsim)):
  dx_PRWPBsim.append(np.diff(np.array(data_PRWPBsim[i]['x'].dropna().tolist())))
  dy_PRWPBsim.append(np.diff(np.array(data_PRWPBsim[i]['y'].dropna().tolist())))

dx_PRWPBsim = np.concatenate(dx_PRWPBsim).ravel()
dy_PRWPBsim = np.concatenate(dy_PRWPBsim).ravel()

dx_dy_PRWPBsim = np.concatenate((dx_PRWPBsim,dy_PRWPBsim))

dx_data = []
dy_data = []
for i in range(len(tracks_geo_region)):
  dx_data.append(np.array(tracks_geo_region[i]['dx'].dropna().tolist()))
  dy_data.append(np.array(tracks_geo_region[i]['dy'].dropna().tolist()))

dx_data = np.concatenate(dx_data).ravel()
dy_data = np.concatenate(dy_data).ravel()

dx_dy_data = np.concatenate((dx_data,dy_data))

# getting data of the histogram
count_data, bins_count_data = np.histogram(dx_dy_data, bins=50)
# finding the PDF of the histogram using count values
pdf_data = count_data / sum(count_data)
# using numpy np.cumsum to calculate the CDF
cdf_data = np.cumsum(pdf_data)

# getting data of the histogram
count_PRWPBsim, bins_count_PRWPBsim = np.histogram(dx_dy_PRWPBsim, bins=50)
# finding the PDF of the histogram using count values
pdf_PRWPBsim = count_PRWPBsim / sum(count_PRWPBsim)
# using numpy np.cumsum to calculate the CDF
cdf_PRWPBsim = np.cumsum(pdf_PRWPBsim)

# getting data of the histogram
count_PRWsim, bins_count_PRWsim = np.histogram(dx_dy_PRWsim, bins=50)
# finding the PDF of the histogram using count values
pdf_PRWsim = count_PRWsim / sum(count_PRWsim)
# using numpy np.cumsum to calculate the CDF
cdf_PRWsim = np.cumsum(pdf_PRWsim)


plt.plot(bins_count_data[1:], cdf_data, label="data")
plt.plot(bins_count_PRWPBsim[1:], cdf_PRWPBsim, label="PRW polarity bias sim, error={}".format(round(np.sum(np.abs(cdf_data-cdf_PRWPBsim)),3)))
plt.plot(bins_count_PRWsim[1:], cdf_PRWsim, label="PRW sim, error={}".format(round(np.sum(np.abs(cdf_data-cdf_PRWsim)),3)))
plt.title('CDF for dx and dy {}'.format(region))
plt.legend()
plt.savefig('figures/model/dxdy_cdf_modelcomaparion_{}'.format(region))
plt.clf()