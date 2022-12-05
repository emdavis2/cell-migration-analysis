from functions.acf_functions import *
from functions.msd_functions import *
from functions.compile_data_tracks_function import *
from functions.PRW_model_functions import *
from functions.langevin_PRW_functions import *
from functions.PRWpolaritybias_model_functions import *
from functions.model_fitting_functions import *

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#from scipy.stats import ttest_ind
from scipy.stats import f_oneway
from scipy.optimize import curve_fit

treatment = str(sys.argv[1])

min_track_length = int(sys.argv[2])

region = str(sys.argv[3])

time = int(sys.argv[4]) #5

dt = float(sys.argv[5]) #0.1667

Nwalkers = int(sys.argv[6]) #113

PRW_params = str(sys.argv[7])

PRWPB_params = str(sys.argv[8])

file_path = str(sys.argv[9])

#get params from model fitting
PRW_params_open = open(PRW_params, 'r')
PRW_params_readlines = PRW_params_open.readlines()
PRW_err = float(PRW_params_readlines[0].split('=')[1])
prw_S = float(PRW_params_readlines[1].split('=')[1])
prw_P = float(PRW_params_readlines[2].split('=')[1])

PRWPB_params_open = open(PRWPB_params, 'r')
PRWPB_params_readlines = PRWPB_params_open.readlines()
PRWPB_err = float(PRWPB_params_readlines[0].split('=')[1])
PRWPB_w_std_dev = float(PRWPB_params_readlines[1].split('=')[1])
PRWPB_theta_std_dev = float(PRWPB_params_readlines[2].split('=')[1])

tracks_region, tracks_geo_region, region_cells, region_endpointcells = compile_data_tracks(treatment, min_track_length, region)

if region == 'stiff':
  region_name = 'gel'
else:
  region_name = region

sampling_t = 10 #min per frame

#clears out sentinel file if it exists
open('sentinels/figs_{}_{}.txt'.format(ntpath.basename(file_path), region),'w').close()
#create new sentinel file to write to
model_data_figs = open('sentinels/figs_{}_{}.txt'.format(file_path[8:-2], region),'w')
file_lines = []

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
data_PRWsim = run_PRW_langevin_sim(Nwalkers, dt, time, prw_S, prw_P)

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
plt.hlines(y=0,xmin=0,xmax=min_track_length,color='k')
plt.xlim(0,min_track_length-2)
plt.ylim(-1,1)
plt.xlabel("time lag")
plt.title(" Autocorrelation velocity {}".format(region_name))
plt.legend()
plt.savefig(file_path + 'acf_velocity_modelcomaparion_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'acf_velocity_modelcomaparion_{}.png \n'.format(region))

#single exponential fits for acf velocity
# Fit the function 
popt, pcov = curve_fit(lambda t, b: np.exp(-t*b), np.arange(0,min_track_length-4),poslagaverage_PRWPBsim[0:min_track_length-4],method='trf')
# Create the fitted curve
b = popt[0]
x_fitted = np.linspace(0, min_track_length-4, 100)
y_fitted = np.exp(-x_fitted*b)
# Plot
plt.scatter(np.arange(0,min_track_length-4), poslagaverage_PRWPBsim[0:min_track_length-4], label='Raw data')
plt.plot(x_fitted, y_fitted, 'k', label='Fitted curve')
plt.title(r'Single exponential fit for acf velocity PRW Polarity Bias {}'.format(region_name))
plt.xlabel(r'time lag $\tau$ (10 min)')
plt.ylabel(r'acf($\tau$)')
plt.legend()
plt.text(12,0.9,r'$acf(\tau)=e^{-\tau b}$')
plt.text(12,0.7,'$b$={}'.format(round(b,3)))
plt.hlines(y=0,xmin=0,xmax=min_track_length,color='k')
plt.xlim(0,min_track_length-2)
plt.ylim(-1,1)
plt.savefig(file_path + 'single_exp_acf_vel_fit_PRWPB_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'single_exp_acf_vel_fit_PRWPB_{}.png \n'.format(region))

# Fit the function 
popt, pcov = curve_fit(lambda t, b: np.exp(-t*b), np.arange(0,min_track_length-4),poslagaverage_data[0:min_track_length-4],method='trf')
# Create the fitted curve
b = popt[0]
x_fitted = np.linspace(0, min_track_length-4, 100)
y_fitted = np.exp(-x_fitted*b)
# Plot
plt.scatter(np.arange(0,min_track_length-4), poslagaverage_data[0:min_track_length-4], label='Raw data')
plt.plot(x_fitted, y_fitted, 'k', label='Fitted curve')
plt.title(r'Single exponential fit for acf velocity {} data'.format(region_name))
plt.xlabel(r'time lag $\tau$ (10 min)')
plt.ylabel(r'acf($\tau$)')
plt.legend()
plt.text(12,0.9,r'$acf(\tau)=e^{-\tau b}$')
plt.text(12,0.7,'$b$={}'.format(round(b,3)))
plt.hlines(y=0,xmin=0,xmax=min_track_length,color='k')
plt.xlim(0,min_track_length-2)
plt.ylim(-1,1)
plt.savefig(file_path + 'single_exp_acf_vel_fit_data_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'single_exp_acf_vel_fit_data_{}.png \n'.format(region))

#double exponential fits for acf velocity
# Fit the function 
popt, pcov = curve_fit(lambda t, b, a: np.exp(-t*b) + np.exp(-t*a), np.arange(0,min_track_length-4),poslagaverage_PRWPBsim[0:min_track_length-4],method='trf')
# Create the fitted curve
b = popt[0]
a = popt[1]
x_fitted = np.linspace(0, min_track_length-4, 100)
y_fitted = np.exp(-x_fitted*b) + np.exp(-x_fitted*a)
# Plot
plt.scatter(np.arange(0,min_track_length-4), poslagaverage_PRWPBsim[0:min_track_length-4], label='Raw data')
plt.plot(x_fitted, y_fitted, 'k', label='Fitted curve')
plt.title(r'Double exponential fit for acf velocity PRW Polarity Bias {}'.format(region_name))
plt.xlabel(r'time lag $\tau$ (10 min)')
plt.ylabel(r'acf($\tau$)')
plt.legend()
plt.text(12,0.9,r'$acf(\tau)=e^{-\tau b} + e^{-\tau b}$')
plt.text(12,0.7,'$b$={} $a$={}'.format(round(b,3),round(a,3)))
plt.hlines(y=0,xmin=0,xmax=min_track_length,color='k')
plt.xlim(0,min_track_length-2)
plt.ylim(-1,1)
plt.savefig(file_path + 'double_exp_acf_vel_fit_PRWPB_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'double_exp_acf_vel_fit_PRWPB_{}.png \n'.format(region))

# Fit the function 
popt, pcov = curve_fit(lambda t, b, a: np.exp(-t*b) + np.exp(-t*a), np.arange(0,min_track_length-4),poslagaverage_data[0:min_track_length-4],method='trf')
# Create the fitted curve
b = popt[0]
a = popt[1]
x_fitted = np.linspace(0, min_track_length-4, 100)
y_fitted = np.exp(-x_fitted*b) + np.exp(-x_fitted*a)
# Plot
plt.scatter(np.arange(0,min_track_length-4), poslagaverage_data[0:min_track_length-4], label='Raw data')
plt.plot(x_fitted, y_fitted, 'k', label='Fitted curve')
plt.title(r'Double exponential fit for acf velocity data {}'.format(region_name))
plt.xlabel(r'time lag $\tau$ (10 min)')
plt.ylabel(r'acf($\tau$)')
plt.legend()
plt.text(12,0.9,r'$acf(\tau)=e^{-\tau b} + e^{-\tau b}$')
plt.text(12,0.7,'$b$={} $a$={}'.format(round(b,3),round(a,3)))
plt.hlines(y=0,xmin=0,xmax=min_track_length,color='k')
plt.xlim(0,min_track_length-2)
plt.ylim(-1,1)
plt.savefig(file_path + 'double_exp_acf_vel_fit_data_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'double_exp_acf_vel_fit_data_{}.png \n'.format(region))

#autocorrelation polarity angle for data
poslagaverage_data = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  normvect = pd.DataFrame( {'cospolangle':list(np.cos(track[['polarity_angle']].iloc[:,0])) , 'sinpolangle':list(np.sin(track[['polarity_angle']].iloc[:,0])) })
  combined = pd.concat([normvect.reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage_data[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
poslagaverage_data /= len(tracks_geo_region) #Nposlagtotal 
poslagaverage_data = poslagaverage_data[0:min_track_length-4]
'''
#autocorrelation polarity angle for PRW model 
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
'''
#autocorrelation polarity angle for PRW polarity bias model
data_PRWPBsim = run_PRWpolaritybias_sim(Nwalkers, dt, time, PRWPB_w_std_dev, PRWPB_theta_std_dev)

poslagaverage_PRWPBsim = np.zeros(300)
all_ac = []
for df in data_PRWPBsim:
    track=df
    combined = make_comb_df(np.cos(track['omega'].to_list()[2:min_track_length-2]),np.sin(track['omega'].to_list()[2:min_track_length-2]))
    combined = combined.dropna()
    poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

    #remove nans here
    poslagsmean[np.isnan(poslagsmean)] = 0
    all_ac.append(poslagsmean)
    poslagaverage_PRWPBsim[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
poslagaverage_PRWPBsim /= len(data_PRWPBsim) #Nposlagtotal 

poslagaverage_PRWPBsim = poslagaverage_PRWPBsim[0:min_track_length-4]

#Compare acf polarity angle plots
plt.plot(poslagaverage_data,label = "Data")
plt.plot(poslagaverage_PRWPBsim,label = "PRW Polarity Bias, error={}".format(round(PRWPB_err,3)))
#plt.plot(poslagaverage_PRWsim,label = "PRW, error={}".format(round(PRW_err,3)))
plt.hlines(y=0,xmin=0,xmax=min_track_length,color='k')
plt.xlim(0,min_track_length-2)
plt.ylim(-1,1)
plt.xlabel("time lag")
plt.title(" Autocorrelation polarity angle {}".format(region_name))
plt.legend()
plt.savefig(file_path + 'acf_polang_modelcomaparion_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'acf_polang_modelcomaparion_{}.png \n'.format(region))

#compare dx and dy with cdf
dx_PRWsim = []
dy_PRWsim = []
for i in range(len(data_PRWsim)):
  dx_PRWsim.append(np.diff(np.array(data_PRWsim[i]['x'].dropna().tolist())))
  dy_PRWsim.append(np.diff(np.array(data_PRWsim[i]['y'].dropna().tolist())))

dx_PRWsim = (np.concatenate(dx_PRWsim).ravel())/sampling_t
dy_PRWsim = (np.concatenate(dy_PRWsim).ravel())/sampling_t

dx_dy_PRWsim = np.concatenate((dx_PRWsim,dy_PRWsim))

dx_PRWPBsim = []
dy_PRWPBsim = []
for i in range(len(data_PRWPBsim)):
  dx_PRWPBsim.append(np.diff(np.array(data_PRWPBsim[i]['x'].dropna().tolist())))
  dy_PRWPBsim.append(np.diff(np.array(data_PRWPBsim[i]['y'].dropna().tolist())))

dx_PRWPBsim = (np.concatenate(dx_PRWPBsim).ravel())/sampling_t
dy_PRWPBsim = (np.concatenate(dy_PRWPBsim).ravel())/sampling_t

dx_dy_PRWPBsim = np.concatenate((dx_PRWPBsim,dy_PRWPBsim))

dx_data = []
dy_data = []
for i in range(len(tracks_geo_region)):
  dx_data.append(np.array(tracks_geo_region[i]['dx'].dropna().tolist()))
  dy_data.append(np.array(tracks_geo_region[i]['dy'].dropna().tolist()))

dx_data = (np.concatenate(dx_data).ravel())/sampling_t
dy_data = (np.concatenate(dy_data).ravel())/sampling_t

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

plt.hist(dx_dy_PRWPBsim,bins=30)
plt.xlabel(r'stepsize ($\mu m$/min)')
plt.title('dx dy for PRW polarity bias sim {}'.format(region_name))
plt.savefig(file_path + 'dxdy_hist_PRWPBsim_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'dxdy_hist_PRWPBsim_{}.png \n'.format(region))

plt.hist(dx_dy_PRWsim,bins=30)
plt.xlabel(r'stepsize ($\mu m$/min)')
plt.title('dx dy for PRW sim {}'.format(region_name))
plt.savefig(file_path + 'dxdy_hist_PRWsim_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'dxdy_hist_PRWsim_{}.png \n'.format(region))

plt.plot(bins_count_data[1:], cdf_data, label="data")
plt.plot(bins_count_PRWPBsim[1:], cdf_PRWPBsim, label="PRW polarity bias sim, error={}".format(round(np.sum(np.abs(cdf_data-cdf_PRWPBsim)),3)))
plt.plot(bins_count_PRWsim[1:], cdf_PRWsim, label="PRW sim, error={}".format(round(np.sum(np.abs(cdf_data-cdf_PRWsim)),3)))
plt.title('CDF for dx and dy {}'.format(region_name))
plt.legend()
plt.savefig(file_path + 'dxdy_cdf_modelcomaparion_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'dxdy_cdf_modelcomaparion_{}.png \n'.format(region))

#Compare dx dy with boxplot
data_bp = {'{} Data'.format(region_name):dx_dy_data, 'PRW Polarity Bias':dx_dy_PRWPBsim, 'PRW':dx_dy_PRWsim}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"DX DY ($\mu m$/min)")
#tstat_PRWPB, pval_PRWPB = ttest_ind(dx_dy_data,dx_dy_PRWPBsim)
tstat_PRWPB, pval_PRWPB = f_oneway(dx_dy_data,dx_dy_PRWPBsim)
plt.text(.2, 40, 'statistic PRW PB={}, pvalue PRW PB={}'.format(round(tstat_PRWPB,3),round(pval_PRWPB,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
#tstat_PRW, pval_PRW = ttest_ind(dx_dy_data,dx_dy_PRWsim)
tstat_PRW, pval_PRW = f_oneway(dx_dy_data,dx_dy_PRWsim)
plt.text(.2, 30, 'statistic PRW={}, pvalue PRW={}'.format(round(tstat_PRW,3),round(pval_PRW,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig(file_path + 'dx_dy_boxplot_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'dx_dy_boxplot_{}.png \n'.format(region))

#Compare vx and vy with boxplot
vx_PRWsim = []
vy_PRWsim = []
v_PRWsim = []
for i in range(len(data_PRWsim)):
  vx_PRWsim.append(np.array(data_PRWsim[i]['vx'].dropna().tolist()))
  vy_PRWsim.append(np.array(data_PRWsim[i]['vy'].dropna().tolist()))
  v_PRWsim.append(np.array(data_PRWsim[i]['v'].dropna().tolist()))

vx_PRWsim = (np.concatenate(vx_PRWsim).ravel())/sampling_t
vy_PRWsim = (np.concatenate(vy_PRWsim).ravel())/sampling_t
v_PRWsim = (np.concatenate(v_PRWsim).ravel())/sampling_t

vx_vy_PRWsim = np.concatenate((vx_PRWsim,vy_PRWsim))

vx_PRWPBsim = []
vy_PRWPBsim = []
v_PRWPBsim = []
for i in range(len(data_PRWPBsim)):
  vx_PRWPBsim.append(np.array(data_PRWPBsim[i]['vx'].dropna().tolist()))
  vy_PRWPBsim.append(np.array(data_PRWPBsim[i]['vy'].dropna().tolist()))
  v_PRWPBsim.append(np.array(data_PRWPBsim[i]['v'].dropna().tolist()))

vx_PRWPBsim = (np.concatenate(vx_PRWPBsim).ravel())/sampling_t
vy_PRWPBsim = (np.concatenate(vy_PRWPBsim).ravel())/sampling_t
v_PRWPBsim = (np.concatenate(v_PRWPBsim).ravel())/sampling_t

vx_vy_PRWPBsim = np.concatenate((vx_PRWPBsim,vy_PRWPBsim))

vx_data = []
vy_data = []
v_data = []
for i in range(len(tracks_geo_region)):
  vx_data.append(np.array(tracks_geo_region[i]['vx'].dropna().tolist()))
  vy_data.append(np.array(tracks_geo_region[i]['vy'].dropna().tolist()))
  v_data.append(np.array(tracks_geo_region[i]['v'].dropna().tolist()))

vx_data = (np.concatenate(vx_data).ravel())/sampling_t
vy_data = (np.concatenate(vy_data).ravel())/sampling_t
v_data = (np.concatenate(v_data).ravel())/sampling_t

vx_vy_data = np.concatenate((vx_data,vy_data))

data_bp = {'{} Data'.format(region_name):vx_vy_data, 'PRW Polarity Bias':vx_vy_PRWPBsim, 'PRW':vx_vy_PRWsim}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"VX VY ($\mu m$/min)")
#tstat_PRWPB, pval_PRWPB = ttest_ind(vx_vy_data,vx_vy_PRWPBsim)
tstat_PRWPB, pval_PRWPB = f_oneway(vx_vy_data,vx_vy_PRWPBsim)
plt.text(.2, -16, 'statistic PRW PB={}, pvalue PRW PB={}'.format(round(tstat_PRWPB,3),round(pval_PRWPB,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
#tstat_PRW, pval_PRW = ttest_ind(vx_vy_data,vx_vy_PRWsim)
tstat_PRW, pval_PRW = f_oneway(vx_vy_data,vx_vy_PRWsim)
plt.text(.2, -20, 'statistic PRW={}, pvalue PRW={}'.format(round(tstat_PRW,3),round(pval_PRW,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig(file_path + 'vx_vy_boxplot_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'vx_vy_boxplot_{}.png \n'.format(region))

data_bp = {'{} Data'.format(region_name):v_data, 'PRW Polarity Bias':v_PRWPBsim, 'PRW':v_PRWsim}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Velocity ($\mu m$/min)")
#tstat_PRWPB, pval_PRWPB = ttest_ind(v_data,v_PRWPBsim)
tstat_PRWPB, pval_PRWPB = f_oneway(v_data,v_PRWPBsim)
plt.text(.2, 22, 'statistic PRW PB={}, pvalue PRW PB={}'.format(round(tstat_PRWPB,3),round(pval_PRWPB,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
#tstat_PRW, pval_PRW = ttest_ind(v_data,v_PRWsim)
tstat_PRW, pval_PRW = f_oneway(v_data,v_PRWsim)
plt.text(.2, 20, 'statistic PRW={}, pvalue PRW={}'.format(round(tstat_PRW,3),round(pval_PRW,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig(file_path + 'velocity_boxplot_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'velocity_boxplot_{}.png \n'.format(region))


#Compare D/T with boxplot
DoverT_PRWsim = []
for i in range(len(data_PRWsim)):
  DoverT_PRWsim.append(np.array(data_PRWsim[i]['DoverT'].dropna().tolist()))
DoverT_PRWsim = np.concatenate(DoverT_PRWsim).ravel()

DoverT_PRWPBsim = []
for i in range(len(data_PRWPBsim)):
  DoverT_PRWPBsim.append(np.array(data_PRWPBsim[i]['DoverT'].dropna().tolist()))
DoverT_PRWPBsim = np.concatenate(DoverT_PRWPBsim).ravel()

DoverT_data = region_endpointcells['DoverT']

data_bp = {'{} Data'.format(region_name):DoverT_data, 'PRW Polarity Bias':DoverT_PRWPBsim, 'PRW':DoverT_PRWsim}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("D/T")
#tstat_PRWPB, pval_PRWPB = ttest_ind(DoverT_data,DoverT_PRWPBsim)
tstat_PRWPB, pval_PRWPB = f_oneway(DoverT_data,DoverT_PRWPBsim)
plt.text(.2, 0.8, 'statistic PRW PB={}, pvalue PRW PB={}'.format(round(tstat_PRWPB,3),round(pval_PRWPB,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
#tstat_PRW, pval_PRW = ttest_ind(DoverT_data,DoverT_PRWsim)
tstat_PRW, pval_PRW = f_oneway(DoverT_data,DoverT_PRWsim)
plt.text(.2, 0.7, 'statistic PRW={}, pvalue PRW={}'.format(round(tstat_PRW,3),round(pval_PRW,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig(file_path + 'DoverT_boxplot_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'DoverT_boxplot_{}.png \n'.format(region))

#Plot trajectories 
#x_data = {}
#y_data = {}
#i = 1
for df in tracks_region:
    plt.plot(df[:,0],df[:,1])
    #x_data[i] = df[:,0]
    #y_data[i] = df[:,1]
    #i += 1
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for {} Data'.format(region_name))
plt.savefig(file_path + 'trajectories_data_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'trajectories_data_{}.png \n'.format(region))

for df in data_PRWPBsim:
  plt.plot(df['x'],df['y'])
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for PRW Polarity Bias')
plt.savefig(file_path + 'trajectories_PRWPBsim_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'trajectories_PRWPBsim_{}.png \n'.format(region))

for df in data_PRWsim:
  plt.plot(df['x'],df['y'])
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for PRW')
plt.savefig(file_path + 'trajectories_PRWsim_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'trajectories_PRWsim_{}.png \n'.format(region))

#track lengths of data
lengths = [len(tracks_geo_region[i]) for i in range(len(tracks_geo_region)) ]
plt.hist(lengths,bins=30)
plt.xlabel('track length (number of frames)')
plt.ylabel('counts')
plt.title('Track lengths for {} data'.format(region_name))
plt.savefig(file_path + 'track_lengths_data_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'track_lengths_data_{}.png \n'.format(region))

#Plotting angles over time
for i in data_PRWPBsim:
  plt.plot(i['omega']-i['omega'][0])
plt.xlabel('time (10 min)')
plt.ylabel('angle (radians)')
plt.title(r'$\omega$ vs time for PRW Polarity Bias Sim {}'.format(region_name))
plt.savefig(file_path + 'omega_vs_time_PRWPB_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'omega_vs_time_PRWPB_{}.png \n'.format(region))

for i in data_PRWPBsim:
  plt.plot(i['theta']-i['theta'][0])
plt.xlabel('time (10 min)')
plt.ylabel('angle (radians)')
plt.title(r'$\theta$ vs time for PRW Polarity Bias Sim {}'.format(region_name))
plt.savefig(file_path + 'theta_vs_time_PRWPB_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'theta_vs_time_PRWPB_{}.png \n'.format(region))

for df in tracks_geo_region:
    plt.plot(df['theta']-df['theta'][0])
plt.xlim(0,min_track_length)
plt.xlabel('time (10 min)')
plt.ylabel('angle (radians)')
plt.title(r'$\theta$ vs time for Data {}'.format(region_name))
plt.savefig(file_path + 'theta_vs_time_data_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'theta_vs_time_data_{}.png \n'.format(region))

#Plot MSD
plt.plot(calc_MSD(tracks_region, min_track_length), label='{} Data'.format(region_name))
plt.plot(calc_MSD_sim(data_PRWPBsim, min_track_length), label='PRW Polarity Bias Sim')
plt.plot(calc_MSD_sim(data_PRWsim, min_track_length), label='PRW Sim')
plt.xlabel('lag (10 min)')
plt.ylabel(r'MSD ($\mu m^{2}$)')
plt.legend()
plt.savefig(file_path + 'MSD_{}.png'.format(region))
plt.clf()
file_lines.append(file_path + 'MSD_{}.png \n'.format(region))

#write lines to text file 
model_data_figs.writelines(file_lines)
model_data_figs.close() 