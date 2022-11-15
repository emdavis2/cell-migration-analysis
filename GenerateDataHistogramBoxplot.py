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

min_track_length = int(sys.argv[3])

region1 = str(sys.argv[4])

region2 = str(sys.argv[5])

tracks_region1, tracks_geo_region1, region1_cells, region1_endpointcells = compile_data_tracks(treatment1, min_track_length, region1)

tracks_region2, tracks_geo_region2, region2_cells, region2_endpointcells = compile_data_tracks(treatment2, min_track_length, region2)

if region2 == 'stiff':
  region2_name = 'gel'
  region1_name = region1
elif region1 == 'stiff':
  region1_name = 'gel'
  region2_name = region2

sampling_t = 10 #min per frame

#Plot histogram for D/T
plt.hist(region1_endpointcells['DoverT'],bins=30)
plt.title('D/T for {}'.format(region1_name))
plt.xlabel('D/T')
plt.ylabel('counts')
plt.savefig('figures/histogram_boxplot/DT_hist_{}.png'.format(region1))
plt.clf()

plt.hist(region2_endpointcells['DoverT'],bins=30)
plt.title('D/T for {}'.format(region2_name))
plt.xlabel('D/T')
plt.ylabel('counts')
plt.savefig('figures/histogram_boxplot/DT_hist_{}.png'.format(region2))
plt.clf()

#Plot histogram for speed
plt.hist(region1_endpointcells['speed']/sampling_t,bins=30)
plt.title('speed for {}'.format(region1_name))
plt.xlabel(r'speed ($\mu m$)')
plt.ylabel('counts')
plt.savefig('figures/histogram_boxplot/speed_hist_{}.png'.format(region1))
plt.clf()

plt.hist(region2_endpointcells['speed']/sampling_t,bins=30)
plt.title('speed for {}'.format(region2_name))
plt.xlabel(r'speed ($\mu m$)')
plt.ylabel('counts')
plt.savefig('figures/histogram_boxplot/speed_hist_{}.png'.format(region2))
plt.clf()

#Plot histogram for FMI
plt.hist(region1_endpointcells['FMI'],bins=30)
plt.title('FMI for {}'.format(region1_name))
plt.xlabel('FMI')
plt.ylabel('counts')
plt.savefig('figures/histogram_boxplot/FMI_hist_{}.png'.format(region1))
plt.clf()

plt.hist(region2_endpointcells['FMI'],bins=30)
plt.title('FMI for {}'.format(region2_name))
plt.xlabel('FMI')
plt.ylabel('counts')
plt.savefig('figures/histogram_boxplot/FMI_hist_{}.png'.format(region2))
plt.clf()


#Plot histogram for velocity
v_region1 = []
for i in range(len(tracks_geo_region1)):
  v_region1.append(tracks_geo_region1[i]['v'].dropna().tolist())

v_region1 = (np.concatenate(v_region1).ravel())/sampling_t

plt.hist(v_region1,bins=50)
plt.title(r'velocity ($\mu m$) {}'.format(region1_name))
plt.savefig('figures/histogram_boxplot/vel_hist_{}.png'.format(region1))
plt.clf()

v_region2 = []
for i in range(len(tracks_geo_region2)):
  v_region2.append(tracks_geo_region2[i]['v'].dropna().tolist())

v_region2 = (np.concatenate(v_region2).ravel())/sampling_t

plt.hist(v_region2,bins=50)
plt.title(r'velocity ($\mu m$) {}'.format(region2_name))
plt.savefig('figures/histogram_boxplot/vel_hist_{}.png'.format(region2))
plt.clf()

#Plot histogram for vx and vy
vx_region1 = []
vy_region1 = []
for i in range(len(tracks_geo_region1)):
  vx_region1.append(tracks_geo_region1[i]['vx'].dropna().tolist())
  vy_region1.append(tracks_geo_region1[i]['vy'].dropna().tolist())

vx_region1 = (np.concatenate(vx_region1).ravel())/sampling_t
vy_region1 = (np.concatenate(vy_region1).ravel())/sampling_t

plt.hist(vx_region1,bins=50)
plt.title(r'vx ($\mu m$) {}'.format(region1_name))
plt.savefig('figures/histogram_boxplot/velx_hist_{}.png'.format(region1))
plt.clf()

plt.hist(vy_region1,bins=50)
plt.title(r'vy ($\mu m$) {}'.format(region1_name))
plt.savefig('figures/histogram_boxplot/vely_hist_{}.png'.format(region1))
plt.clf()

vx_region2 = []
vy_region2 = []
for i in range(len(tracks_geo_region2)):
  vx_region2.append(tracks_geo_region2[i]['vx'].dropna().tolist())
  vy_region2.append(tracks_geo_region2[i]['vy'].dropna().tolist())

vx_region2 = (np.concatenate(vx_region2).ravel())/sampling_t
vy_region2 = (np.concatenate(vy_region2).ravel())/sampling_t

plt.hist(vx_region2,bins=50)
plt.title(r'vx ($\mu m$) {}'.format(region2_name))
plt.savefig('figures/histogram_boxplot/velx_hist_{}.png'.format(region2))
plt.clf()

plt.hist(vy_region2,bins=50)
plt.title(r'vy ($\mu m$) {}'.format(region2_name))
plt.savefig('figures/histogram_boxplot/vely_hist_{}.png'.format(region2))
plt.clf()

#Plot histogram for abs-skew
absskew_region1 = []
for i in range(len(tracks_geo_region1)):
  absskew_region1.append(tracks_geo_region1[i]['abs-skew'].dropna().tolist())

absskew_region1 = np.concatenate(absskew_region1).ravel()

plt.hist(absskew_region1,bins=50)
plt.title('abs-skew {}'.format(region1_name))
plt.savefig('figures/histogram_boxplot/absskew_hist_{}.png'.format(region1))
plt.clf()

absskew_region2 = []
for i in range(len(tracks_geo_region2)):
  absskew_region2.append(tracks_geo_region2[i]['abs-skew'].dropna().tolist())

absskew_region2 = np.concatenate(absskew_region2).ravel()

plt.hist(absskew_region2,bins=50)
plt.title('abs-skew {}'.format(region2_name))
plt.savefig('figures/histogram_boxplot/absskew_hist_{}.png'.format(region2))
plt.clf()

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
plt.title(r'dx dy ($\mu m$) {}'.format(region1_name))
plt.savefig('figures/histogram_boxplot/dxdy_hist_{}.png'.format(region1))
plt.clf()

plt.hist(dx_region1,bins=50)
plt.title(r'dx ($\mu m$) {}'.format(region1_name))
plt.savefig('figures/histogram_boxplot/dx_hist_{}.png'.format(region1))
plt.clf()

plt.hist(dy_region1,bins=50)
plt.title(r'dy ($\mu m$) {}'.format(region1_name))
plt.savefig('figures/histogram_boxplot/dy_hist_{}.png'.format(region1))
plt.clf()

dx_region2 = []
dy_region2 = []
for i in range(len(tracks_geo_region2)):
  dx_region2.append(tracks_geo_region2[i]['dx'].dropna())
  dy_region2.append(tracks_geo_region2[i]['dy'].dropna())

dx_region2 = (np.concatenate(dx_region2).ravel())/sampling_t
dy_region2 = (np.concatenate(dy_region2).ravel())/sampling_t

dx_dy_region2 = np.concatenate((dx_region2,dy_region2))

plt.hist(dx_dy_region2,bins=50)
plt.title(r'dx dy ($\mu m$) {}'.format(region2_name))
plt.savefig('figures/histogram_boxplot/dxdy_hist_{}.png'.format(region2))
plt.clf()

plt.hist(dx_region2,bins=50)
plt.title(r'dx ($\mu m$) {}'.format(region2_name))
plt.savefig('figures/histogram_boxplot/dx_hist_{}.png'.format(region2))
plt.clf()

plt.hist(dy_region2,bins=50)
plt.title(r'dy ($\mu m$) {}'.format(region2_name))
plt.savefig('figures/histogram_boxplot/dy_hist_{}.png'.format(region2))
plt.clf()

#make boxplots
data_bp = {'{}'.format(region1_name):v_region1, '{}'.format(region2_name):v_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Velocity ($\mu m$)")
#tstat, pval = ttest_ind(v_region1,v_region2)
tstat, pval = f_oneway(v_region1,v_region2)
plt.text(.2, 20, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/velocity_boxplot.png')
plt.clf()

data_bp = {'{}'.format(region1_name):vx_region1, '{}'.format(region2_name):vx_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Velocity x ($\mu m$)")
#tstat, pval = ttest_ind(vx_region1,vx_region2)
tstat, pval = f_oneway(vx_region1,vx_region2)
plt.text(.1, 18, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/velocity_x_boxplot.png')
plt.clf()

data_bp = {'{}'.format(region1_name):vy_region1, '{}'.format(region2_name):vy_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Velocity y ($\mu m$)")
#tstat, pval = ttest_ind(vy_region1,vy_region2)
tstat, pval = f_oneway(vy_region1,vy_region2)
plt.text(.1, 12, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/velocity_y_boxplot.png')
plt.clf()

data_bp = {'{}'.format(region1_name):absskew_region1, '{}'.format(region2_name):absskew_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("abs-skew")
#tstat, pval = ttest_ind(absskew_region1,absskew_region2)
tstat, pval = f_oneway(absskew_region1,absskew_region2)
plt.text(2, 12, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/absskew_boxplot.png')
plt.clf()

data_bp = {'{}'.format(region1_name):dx_region1, '{}'.format(region2_name):dx_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"dx ($\mu m$)")
#tstat, pval = ttest_ind(dx_region1,dx_region2)
tstat, pval = f_oneway(dx_region1,dx_region2)
plt.text(.1, 12, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/dx_boxplot.png')
plt.clf()

data_bp = {'{}'.format(region1_name):dy_region1, '{}'.format(region2_name):dy_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"dy ($\mu m$)")
#tstat, pval = ttest_ind(dy_region1,dy_region2)
tstat, pval = f_oneway(dy_region1,dy_region2)
plt.text(.1, 12, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/dy_boxplot.png')
plt.clf()

data_bp = {'{}'.format(region1_name):region1_endpointcells['speed'], '{}'.format(region2_name):region2_endpointcells['speed']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Speed ($\mu m$)")
#tstat, pval = ttest_ind(region1_endpointcells['speed'],region2_endpointcells['speed'])
tstat, pval = f_oneway(region1_endpointcells['speed'],region2_endpointcells['speed'])
plt.text(.1, 15, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/speed_boxplot.png')
plt.clf()

data_bp = {'{}'.format(region1_name):region1_endpointcells['DoverT'], '{}'.format(region2_name):region2_endpointcells['DoverT']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("D/T")
#tstat, pval = ttest_ind(region1_endpointcells['DoverT'],region2_endpointcells['DoverT'])
tstat, pval = f_oneway(region1_endpointcells['DoverT'],region2_endpointcells['DoverT'])
plt.text(.1, 0.8, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/DoverT_boxplot.png')
plt.clf()

data_bp = {'{}'.format(region1_name):region1_endpointcells['FMI'], '{}'.format(region2_name):region2_endpointcells['FMI']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("FMI")
#tstat, pval = ttest_ind(region1_endpointcells['FMI'],region2_endpointcells['FMI'])
tstat, pval = f_oneway(region1_endpointcells['FMI'],region2_endpointcells['FMI'])
plt.text(.1, 0.7, 'statistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/histogram_boxplot/FMI_boxplot.png')
plt.clf()