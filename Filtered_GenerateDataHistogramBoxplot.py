
from functions.compile_data_tracks_function import *

import ntpath
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
#from scipy.stats import ttest_ind
from scipy.stats import f_oneway

treatment1 = str(sys.argv[1])

treatment2 = str(sys.argv[2])

treatment3 = str(sys.argv[3])

min_track_length = int(sys.argv[4])

region1 = str(sys.argv[5])

region2 = str(sys.argv[6])

region3 = str(sys.argv[7])

pixel_size = 1.54

min_persistence = 0.6

tracks_region1, tracks_geo_region1, region1_cells, region1_endpointcells = compile_data_tracks(treatment1, min_track_length, region1, pixel_size)

tracks_region2, tracks_geo_region2, region2_cells, region2_endpointcells = compile_data_tracks(treatment2, min_track_length, region2, pixel_size)

tracks_region3, tracks_geo_region3, region3_cells, region3_endpointcells = compile_data_tracks(treatment3, min_track_length, region3, pixel_size)

if region2 == 'stiff' or region2 == 'gel':
  region2_name = 'gel'
  region1_name = region1
  region3_name = region3
elif region1 == 'stiff' or region1 == 'gel':
  region1_name = 'gel'
  region2_name = region2
  region3_name = region3
elif region3 == 'stiff' or region3 == 'gel':
  region3_name = 'gel'
  region1_name = region1
  region2_name = region2

else:
  region1_name = region1
  region2_name = region2
  region3_name = region3


sampling_t = 5 #min per frame

#clears out sentinel file if it exists
open('sentinels/filtered_histogram_boxplot.txt','w').close()
#create new sentinel file to write to
hist_boxplot_figs = open('sentinels/filtered_histogram_boxplot.txt','w')
file_lines = []

#filter based on D/T
region1_endpointcells = region1_endpointcells[region1_endpointcells['DoverT'] > min_persistence] 
region2_endpointcells = region2_endpointcells[region2_endpointcells['DoverT'] > min_persistence] 
region3_endpointcells = region3_endpointcells[region3_endpointcells['DoverT'] > min_persistence] 

#names of experiments that meet this requirement
region1_experiments = []
for i in range(len(region1_endpointcells)):
  key_name = str(region1_endpointcells.iloc[i]['experiment']) + '/' + str(region1_endpointcells.iloc[i]['track_id'])
  region1_experiments.append(key_name)

region2_experiments = []
for i in range(len(region2_endpointcells)):
  key_name = str(region2_endpointcells.iloc[i]['experiment']) + '/' + str(region2_endpointcells.iloc[i]['track_id'])
  region2_experiments.append(key_name)

region3_experiments = []
for i in range(len(region3_endpointcells)):
  key_name = str(region3_endpointcells.iloc[i]['experiment']) + '/' + str(region3_endpointcells.iloc[i]['track_id'])
  region3_experiments.append(key_name)



filtered_tracks_geo_region1 = []
filtered_tracks_region1 = []
for i in range(len(tracks_geo_region1)):
  trial_df = str(tracks_geo_region1[i]['experiment'][0]) + '/' + str(tracks_geo_region1[i]['track_id'][0])
  if trial_df in region1_experiments:
    filtered_tracks_geo_region1.append(tracks_geo_region1[i])
    filtered_tracks_region1.append(tracks_region1[i])
  

filtered_tracks_geo_region2 = []
filtered_tracks_region2 = []
for i in range(len(tracks_geo_region2)):
  trial_df = str(tracks_geo_region2[i]['experiment'][0]) + '/' + str(tracks_geo_region2[i]['track_id'][0])
  if trial_df in region2_experiments:
    filtered_tracks_geo_region2.append(tracks_geo_region2[i])
    filtered_tracks_region2.append(tracks_region2[i])

filtered_tracks_geo_region3 = []
filtered_tracks_region3 = []
for i in range(len(tracks_geo_region3)):
  trial_df = str(tracks_geo_region3[i]['experiment'][0]) + '/' + str(tracks_geo_region3[i]['track_id'][0])
  if trial_df in region3_experiments:
    filtered_tracks_geo_region3.append(tracks_geo_region3[i])
    filtered_tracks_region3.append(tracks_region3[i])

tracks_geo_region1 = filtered_tracks_geo_region1
tracks_geo_region2 = filtered_tracks_geo_region2
tracks_geo_region3 = filtered_tracks_geo_region3

tracks_region1 = filtered_tracks_region1
tracks_region2 = filtered_tracks_region2
tracks_region3 = filtered_tracks_region3

#Plot histogram for D/T
plt.hist(region1_endpointcells['DoverT'],bins=30)
plt.title('D/T for {}'.format(region1_name))
plt.xlabel('D/T')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/DT_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/DT_hist_{}.png \n'.format(region1))

plt.hist(region2_endpointcells['DoverT'],bins=30)
plt.title('D/T for {}'.format(region2_name))
plt.xlabel('D/T')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/DT_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/DT_hist_{}.png \n'.format(region2))

plt.hist(region3_endpointcells['DoverT'],bins=30)
plt.title('D/T for {}'.format(region3_name))
plt.xlabel('D/T')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/DT_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/DT_hist_{}.png \n'.format(region3))

#Plot histogram for speed
plt.hist(region1_endpointcells['speed']/sampling_t,bins=30)
plt.title('speed for {}'.format(region1_name))
plt.xlabel(r'speed ($\mu m$/min)')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/speed_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/speed_hist_{}.png \n'.format(region1))

plt.hist(region2_endpointcells['speed']/sampling_t,bins=30)
plt.title('speed for {}'.format(region2_name))
plt.xlabel(r'speed ($\mu m$/min)')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/speed_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/speed_hist_{}.png \n'.format(region2))

plt.hist(region3_endpointcells['speed']/sampling_t,bins=30)
plt.title('speed for {}'.format(region3_name))
plt.xlabel(r'speed ($\mu m$/min)')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/speed_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/speed_hist_{}.png \n'.format(region3))

#Plot histogram for FMI
plt.hist(region1_endpointcells['FMI'],bins=30)
plt.title('FMI for {}'.format(region1_name))
plt.xlabel('FMI')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/FMI_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/FMI_hist_{}.png \n'.format(region1))

plt.hist(region2_endpointcells['FMI'],bins=30)
plt.title('FMI for {}'.format(region2_name))
plt.xlabel('FMI')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/FMI_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/FMI_hist_{}.png \n'.format(region2))

plt.hist(region3_endpointcells['FMI'],bins=30)
plt.title('FMI for {}'.format(region3_name))
plt.xlabel('FMI')
plt.ylabel('counts')
plt.savefig('figures/filtered_histogram_boxplot/FMI_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/FMI_hist_{}.png \n'.format(region3))

#Plot histogram for velocity
v_region1 = []
for i in range(len(tracks_geo_region1)):
  v_region1.append(tracks_geo_region1[i]['v'].dropna().tolist())

v_region1 = (np.concatenate(v_region1).ravel())/sampling_t

plt.hist(v_region1,bins=50)
plt.title(r'velocity ($\mu m$/min) {}'.format(region1_name))
plt.savefig('figures/filtered_histogram_boxplot/vel_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/vel_hist_{}.png \n'.format(region1))

v_region2 = []
for i in range(len(tracks_geo_region2)):
  v_region2.append(tracks_geo_region2[i]['v'].dropna().tolist())

v_region2 = (np.concatenate(v_region2).ravel())/sampling_t

plt.hist(v_region2,bins=50)
plt.title(r'velocity ($\mu m$/min) {}'.format(region2_name))
plt.savefig('figures/filtered_histogram_boxplot/vel_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/vel_hist_{}.png \n'.format(region2))

v_region3 = []
for i in range(len(tracks_geo_region3)):
  v_region3.append(tracks_geo_region3[i]['v'].dropna().tolist())

v_region3 = (np.concatenate(v_region3).ravel())/sampling_t

plt.hist(v_region3,bins=50)
plt.title(r'velocity ($\mu m$/min) {}'.format(region3_name))
plt.savefig('figures/filtered_histogram_boxplot/vel_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/vel_hist_{}.png \n'.format(region3))

#Plot histogram for vx and vy
vx_region1 = []
vy_region1 = []
for i in range(len(tracks_geo_region1)):
  vx_region1.append(tracks_geo_region1[i]['vx'].dropna().tolist())
  vy_region1.append(tracks_geo_region1[i]['vy'].dropna().tolist())

vx_region1 = (np.concatenate(vx_region1).ravel())/sampling_t
vy_region1 = (np.concatenate(vy_region1).ravel())/sampling_t

plt.hist(vx_region1,bins=50)
plt.title(r'vx ($\mu m$/min) {}'.format(region1_name))
plt.savefig('figures/filtered_histogram_boxplot/velx_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/velx_hist_{}.png \n'.format(region1))

plt.hist(vy_region1,bins=50)
plt.title(r'vy ($\mu m$/min) {}'.format(region1_name))
plt.savefig('figures/filtered_histogram_boxplot/vely_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/vely_hist_{}.png \n'.format(region1))

vx_region2 = []
vy_region2 = []
for i in range(len(tracks_geo_region2)):
  vx_region2.append(tracks_geo_region2[i]['vx'].dropna().tolist())
  vy_region2.append(tracks_geo_region2[i]['vy'].dropna().tolist())

vx_region2 = (np.concatenate(vx_region2).ravel())/sampling_t
vy_region2 = (np.concatenate(vy_region2).ravel())/sampling_t

plt.hist(vx_region2,bins=50)
plt.title(r'vx ($\mu m$/min) {}'.format(region2_name))
plt.savefig('figures/filtered_histogram_boxplot/velx_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/velx_hist_{}.png \n'.format(region2))

plt.hist(vy_region2,bins=50)
plt.title(r'vy ($\mu m$/min) {}'.format(region2_name))
plt.savefig('figures/filtered_histogram_boxplot/vely_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/vely_hist_{}.png \n'.format(region2))

vx_region3 = []
vy_region3 = []
for i in range(len(tracks_geo_region3)):
  vx_region3.append(tracks_geo_region3[i]['vx'].dropna().tolist())
  vy_region3.append(tracks_geo_region3[i]['vy'].dropna().tolist())

vx_region3 = (np.concatenate(vx_region3).ravel())/sampling_t
vy_region3 = (np.concatenate(vy_region3).ravel())/sampling_t

plt.hist(vx_region3,bins=50)
plt.title(r'vx ($\mu m$/min) {}'.format(region3_name))
plt.savefig('figures/filtered_histogram_boxplot/velx_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/velx_hist_{}.png \n'.format(region3))

plt.hist(vy_region3,bins=50)
plt.title(r'vy ($\mu m$/min) {}'.format(region3_name))
plt.savefig('figures/filtered_histogram_boxplot/vely_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/vely_hist_{}.png \n'.format(region3))

#Plot histogram for abs-skew
absskew_region1 = []
for i in range(len(tracks_geo_region1)):
  absskew_region1.append(tracks_geo_region1[i]['abs-skew'].dropna().tolist())

absskew_region1 = np.concatenate(absskew_region1).ravel()

plt.hist(absskew_region1,bins=50)
plt.title('abs-skew {}'.format(region1_name))
plt.savefig('figures/filtered_histogram_boxplot/absskew_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/absskew_hist_{}.png \n'.format(region1))

absskew_region2 = []
for i in range(len(tracks_geo_region2)):
  absskew_region2.append(tracks_geo_region2[i]['abs-skew'].dropna().tolist())

absskew_region2 = np.concatenate(absskew_region2).ravel()

plt.hist(absskew_region2,bins=50)
plt.title('abs-skew {}'.format(region2_name))
plt.savefig('figures/filtered_histogram_boxplot/absskew_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/absskew_hist_{}.png \n'.format(region2))

absskew_region3 = []
for i in range(len(tracks_geo_region3)):
  absskew_region3.append(tracks_geo_region3[i]['abs-skew'].dropna().tolist())

absskew_region3 = np.concatenate(absskew_region3).ravel()

plt.hist(absskew_region3,bins=50)
plt.title('abs-skew {}'.format(region3_name))
plt.savefig('figures/filtered_histogram_boxplot/absskew_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/absskew_hist_{}.png \n'.format(region3))

#Plot histogram for dx and dy
dx_region1 = []
dy_region1 = []
for i in range(len(tracks_geo_region1)):
  dx_region1.append(tracks_geo_region1[i]['dx'].dropna())
  dy_region1.append(tracks_geo_region1[i]['dy'].dropna())

dx_region1 = (np.concatenate(dx_region1).ravel())/sampling_t
dy_region1 = (np.concatenate(dy_region1).ravel())/sampling_t

dx_dy_region1 = np.concatenate((dx_region1,dy_region1))

plt.hist(dx_dy_region1,bins=50)
plt.title(r'dx dy ($\mu m$/min) {}'.format(region1_name))
plt.savefig('figures/filtered_histogram_boxplot/dxdy_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dxdy_hist_{}.png \n'.format(region1))

plt.hist(dx_region1,bins=50)
plt.title(r'dx ($\mu m$/min) {}'.format(region1_name))
plt.savefig('figures/filtered_histogram_boxplot/dx_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dx_hist_{}.png \n'.format(region1))

plt.hist(dy_region1,bins=50)
plt.title(r'dy ($\mu m$/min) {}'.format(region1_name))
plt.savefig('figures/filtered_histogram_boxplot/dy_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dy_hist_{}.png \n'.format(region1))

dx_region2 = []
dy_region2 = []
for i in range(len(tracks_geo_region2)):
  dx_region2.append(tracks_geo_region2[i]['dx'].dropna())
  dy_region2.append(tracks_geo_region2[i]['dy'].dropna())

dx_region2 = (np.concatenate(dx_region2).ravel())/sampling_t
dy_region2 = (np.concatenate(dy_region2).ravel())/sampling_t

dx_dy_region2 = np.concatenate((dx_region2,dy_region2))

plt.hist(dx_dy_region2,bins=50)
plt.title(r'dx dy ($\mu m$/min) {}'.format(region2_name))
plt.savefig('figures/filtered_histogram_boxplot/dxdy_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dxdy_hist_{}.png \n'.format(region2))

plt.hist(dx_region2,bins=50)
plt.title(r'dx ($\mu m$/min) {}'.format(region2_name))
plt.savefig('figures/filtered_histogram_boxplot/dx_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dx_hist_{}.png \n'.format(region2))

plt.hist(dy_region2,bins=50)
plt.title(r'dy ($\mu m$/min) {}'.format(region2_name))
plt.savefig('figures/filtered_histogram_boxplot/dy_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dy_hist_{}.png \n'.format(region2))

dx_region3 = []
dy_region3 = []
for i in range(len(tracks_geo_region3)):
  dx_region3.append(tracks_geo_region3[i]['dx'].dropna())
  dy_region3.append(tracks_geo_region3[i]['dy'].dropna())

dx_region3 = (np.concatenate(dx_region3).ravel())/sampling_t
dy_region3 = (np.concatenate(dy_region3).ravel())/sampling_t

dx_dy_region3 = np.concatenate((dx_region3,dy_region3))

plt.hist(dx_dy_region3,bins=50)
plt.title(r'dx dy ($\mu m$/min) {}'.format(region3_name))
plt.savefig('figures/filtered_histogram_boxplot/dxdy_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dxdy_hist_{}.png \n'.format(region3))

plt.hist(dx_region3,bins=50)
plt.title(r'dx ($\mu m$/min) {}'.format(region3_name))
plt.savefig('figures/filtered_histogram_boxplot/dx_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dx_hist_{}.png \n'.format(region3))

plt.hist(dy_region3,bins=50)
plt.title(r'dy ($\mu m$/min) {}'.format(region3_name))
plt.savefig('figures/filtered_histogram_boxplot/dy_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dy_hist_{}.png \n'.format(region3))

#Get solidity
solidity_region1 = []
for i in range(len(tracks_geo_region1)):
  solidity_region1.append(tracks_geo_region1[i]['solidity'].tolist())

solidity_region1 = np.concatenate(solidity_region1).ravel()

solidity_region2 = []
for i in range(len(tracks_geo_region2)):
  solidity_region2.append(tracks_geo_region2[i]['solidity'].tolist())

solidity_region2 = np.concatenate(solidity_region2).ravel()

solidity_region3 = []
for i in range(len(tracks_geo_region3)):
  solidity_region3.append(tracks_geo_region3[i]['solidity'].tolist())

solidity_region3 = np.concatenate(solidity_region3).ravel()

#make boxplots
data_bp = {'{}'.format(region1_name):solidity_region1, '{}'.format(region2_name):solidity_region2, '{}'.format(region3_name):solidity_region3}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("Solidity")
tstat12, pval12 = f_oneway(solidity_region1,solidity_region2)
tstat13, pval13 = f_oneway(solidity_region1,solidity_region3)
tstat23, pval23 = f_oneway(solidity_region2,solidity_region3)
# plt.text(0.1, 0.8, 'pvalue_12={}'.format(pval12), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(0.1, 0.7, 'pvalue_13={}'.format(pval13), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(0.1, 0.6, 'pvalue_23={}'.format(pval23), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/solidity_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/solidity_boxplot.png \n')

data_bp = {'{}'.format(region1_name):v_region1, '{}'.format(region2_name):v_region2, '{}'.format(region3_name):v_region3}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Velocity ($\mu m$/min)")
#tstat, pval = ttest_ind(v_region1,v_region2)
tstat12, pval12 = f_oneway(v_region1,v_region2)
tstat13, pval13 = f_oneway(v_region1,v_region3)
tstat23, pval23 = f_oneway(v_region2,v_region3)
plt.text(.2, 2, 'statistic_12={}, pvalue_12={}'.format(round(tstat12,3),round(pval12,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.2, 1, 'statistic_13={}, pvalue_13={}'.format(round(tstat13,3),round(pval13,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.2, .2, 'statistic_23={}, pvalue_23={}'.format(round(tstat23,3),round(pval23,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/velocity_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/velocity_boxplot.png \n')

data_bp = {'{}'.format(region1_name):vx_region1, '{}'.format(region2_name):vx_region2, '{}'.format(region3_name):vx_region3}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Velocity x ($\mu m$/min)")
#tstat, pval = ttest_ind(vx_region1,vx_region2)
tstat12, pval12 = f_oneway(vx_region1,vx_region2)
tstat13, pval13 = f_oneway(vx_region1,vx_region3)
tstat23, pval23 = f_oneway(vx_region2,vx_region3)
plt.text(.1, 5, 'statistic_12={}, pvalue_12={}'.format(round(tstat12,3),round(pval12,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 4.5, 'statistic_13={}, pvalue_13={}'.format(round(tstat13,3),round(pval13,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 4, 'statistic_23={}, pvalue_23={}'.format(round(tstat23,3),round(pval23,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/velocity_x_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/velocity_x_boxplot.png \n')

data_bp = {'{}'.format(region1_name):vy_region1, '{}'.format(region2_name):vy_region2, '{}'.format(region3_name):vy_region3}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Velocity y ($\mu m$/min)")
#tstat, pval = ttest_ind(vy_region1,vy_region2)
tstat12, pval12 = f_oneway(vy_region1,vy_region2)
tstat13, pval13 = f_oneway(vy_region1,vy_region3)
tstat23, pval23 = f_oneway(vy_region2,vy_region3)
plt.text(.1, 4, 'statistic_12={}, pvalue_12={}'.format(round(tstat12,3),round(pval12,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 3, 'statistic_13={}, pvalue_13={}'.format(round(tstat13,3),round(pval13,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 2, 'statistic_23={}, pvalue_23={}'.format(round(tstat23,3),round(pval23,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/velocity_y_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/velocity_y_boxplot.png \n')

data_bp = {'{}'.format(region1_name):absskew_region1, '{}'.format(region2_name):absskew_region2, '{}'.format(region3_name):absskew_region3}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("abs-skew")
#tstat, pval = ttest_ind(absskew_region1,absskew_region2)
tstat12, pval12 = f_oneway(absskew_region1,absskew_region2)
tstat13, pval13 = f_oneway(absskew_region1,absskew_region3)
tstat23, pval23 = f_oneway(absskew_region2,absskew_region3)
plt.text(.2, 2, 'statistic_12={}, pvalue_12={}'.format(round(tstat12,3),round(pval12,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.2, 1.5, 'statistic_13={}, pvalue_13={}'.format(round(tstat13,3),round(pval13,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.2, 1, 'statistic_23={}, pvalue_23={}'.format(round(tstat23,3),round(pval23,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/absskew_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/absskew_boxplot.png \n')

data_bp = {'{}'.format(region1_name):dx_region1, '{}'.format(region2_name):dx_region2, '{}'.format(region3_name):dx_region3}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"dx ($\mu m$/min)")
#tstat, pval = ttest_ind(dx_region1,dx_region2)
tstat12, pval12 = f_oneway(dx_region1,dx_region2)
tstat13, pval13 = f_oneway(dx_region1,dx_region3)
tstat23, pval23 = f_oneway(dx_region2,dx_region3)
plt.text(.1, 10, 'statistic_12={}, pvalue_12={}'.format(round(tstat12,3),round(pval12,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 9, 'statistic_13={}, pvalue_13={}'.format(round(tstat13,3),round(pval13,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 8, 'statistic_23={}, pvalue_23={}'.format(round(tstat23,3),round(pval23,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/dx_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dx_boxplot.png \n')

data_bp = {'{}'.format(region1_name):dy_region1, '{}'.format(region2_name):dy_region2, '{}'.format(region3_name):dy_region3}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"dy ($\mu m$/min)")
#tstat, pval = ttest_ind(dy_region1,dy_region2)
tstat12, pval12 = f_oneway(dy_region1,dy_region2)
tstat13, pval13 = f_oneway(dy_region1,dy_region3)
tstat13, pval23 = f_oneway(dy_region2,dy_region3)
plt.text(.1, 9, 'statistic_12={}, pvalue_12={}'.format(round(tstat12,3),round(pval12,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 8, 'statistic_13={}, pvalue_13={}'.format(round(tstat13,3),round(pval13,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 7, 'statistic_23={}, pvalue_23={}'.format(round(tstat23,3),round(pval23,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/dy_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/dy_boxplot.png \n')

data_bp = {'{}'.format(region1_name):region1_endpointcells['speed']/sampling_t, '{}'.format(region2_name):region2_endpointcells['speed']/sampling_t, '{}'.format(region3_name):region3_endpointcells['speed']/sampling_t}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Speed ($\mu m$/min)")
#tstat, pval = ttest_ind(region1_endpointcells['speed'],region2_endpointcells['speed'])
tstat12, pval12 = f_oneway(region1_endpointcells['speed'],region2_endpointcells['speed'])
tstat13, pval13 = f_oneway(region1_endpointcells['speed'],region3_endpointcells['speed'])
tstat23, pval23 = f_oneway(region2_endpointcells['speed'],region3_endpointcells['speed'])
# plt.text(.1, 1.5, 'pvalue_12={}'.format(pval12), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(.1, 1.2, 'pvalue_13={}'.format(pval13), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(.1, 1, 'pvalue_23={}'.format(pval23), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/speed_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/speed_boxplot.png \n')

data_bp = {'{}'.format(region1_name):region1_endpointcells['DoverT'], '{}'.format(region2_name):region2_endpointcells['DoverT'], '{}'.format(region3_name):region3_endpointcells['DoverT']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("D/T")
#tstat, pval = ttest_ind(region1_endpointcells['DoverT'],region2_endpointcells['DoverT'])
tstat12, pval12 = f_oneway(region1_endpointcells['DoverT'],region2_endpointcells['DoverT'])
tstat13, pval13 = f_oneway(region1_endpointcells['DoverT'],region3_endpointcells['DoverT'])
tstat23, pval23 = f_oneway(region2_endpointcells['DoverT'],region3_endpointcells['DoverT'])
# plt.text(.1, 0.4, 'pvalue_12={}'.format(pval12), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(.1, 0.3, 'pvalue_13={}'.format(pval13), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(.1, 0.2, 'pvalue_23={}'.format(pval23), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/DoverT_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/DoverT_boxplot.png \n')

data_bp = {'{}'.format(region1_name):region1_endpointcells['FMI'], '{}'.format(region2_name):region2_endpointcells['FMI'], '{}'.format(region3_name):region3_endpointcells['FMI']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("FMI")
#tstat, pval = ttest_ind(region1_endpointcells['FMI'],region2_endpointcells['FMI'])
tstat12, pval12 = f_oneway(region1_endpointcells['FMI'],region2_endpointcells['FMI'])
tstat13, pval13 = f_oneway(region1_endpointcells['FMI'],region3_endpointcells['FMI'])
tstat23, pval23 = f_oneway(region2_endpointcells['FMI'],region3_endpointcells['FMI'])
plt.text(.1, 0.4, 'statistic_12={}, pvalue_12={}'.format(round(tstat12,3),round(pval12,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 0.3, 'statistic_13={}, pvalue_13={}'.format(round(tstat13,3),round(pval13,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.text(.1, 0.2, 'statistic_23={}, pvalue_23={}'.format(round(tstat23,3),round(pval23,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/filtered_histogram_boxplot/FMI_boxplot.png')
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/FMI_boxplot.png \n')


#track lengths of data
lengths = [len(tracks_geo_region1[i]) for i in range(len(tracks_geo_region1)) ]
num_tracks = len(tracks_geo_region1)
plt.hist(lengths)
plt.xlabel('track length (number of frames)')
plt.ylabel('counts')
plt.title('Track lengths for {} data ({} tracks)'.format(region1_name, num_tracks))
plt.savefig('figures/filtered_histogram_boxplot/track_lengths_data_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/track_lengths_data_{}.png \n'.format(region1))

#track lengths of data
lengths = [len(tracks_geo_region2[i]) for i in range(len(tracks_geo_region2)) ]
num_tracks = len(tracks_geo_region2)
plt.hist(lengths)
plt.xlabel('track length (number of frames)')
plt.ylabel('counts')
plt.title('Track lengths for {} data ({} tracks)'.format(region2_name, num_tracks))
plt.savefig('figures/filtered_histogram_boxplot/track_lengths_data_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/track_lengths_data_{}.png \n'.format(region2))

#track lengths of data
lengths = [len(tracks_geo_region3[i]) for i in range(len(tracks_geo_region3)) ]
num_tracks = len(tracks_geo_region3)
plt.hist(lengths)
plt.xlabel('track length (number of frames)')
plt.ylabel('counts')
plt.title('Track lengths for {} data ({} tracks)'.format(region3_name, num_tracks))
plt.savefig('figures/filtered_histogram_boxplot/track_lengths_data_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/track_lengths_data_{}.png \n'.format(region3))

#Trajectories
for df in tracks_region1:
    plt.plot(df[:,0],df[:,1])
    #x_data[i] = df[:,0]
    #y_data[i] = df[:,1]
    #i += 1
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for {} Data'.format(region1_name))
plt.savefig('figures/filtered_histogram_boxplot/trajectories_data_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/trajectories_data_{}.png \n'.format(region1))

for df in tracks_region2:
    plt.plot(df[:,0],df[:,1])
    #x_data[i] = df[:,0]
    #y_data[i] = df[:,1]
    #i += 1
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for {} Data'.format(region2_name))
plt.savefig('figures/filtered_histogram_boxplot/trajectories_data_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/trajectories_data_{}.png \n'.format(region2))

for df in tracks_region3:
    plt.plot(df[:,0],df[:,1])
    #x_data[i] = df[:,0]
    #y_data[i] = df[:,1]
    #i += 1
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for {} Data'.format(region3_name))
plt.savefig('figures/filtered_histogram_boxplot/trajectories_data_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/filtered_histogram_boxplot/trajectories_data_{}.png \n'.format(region3))

#write lines to text file 
hist_boxplot_figs.writelines(file_lines)
hist_boxplot_figs.close() 