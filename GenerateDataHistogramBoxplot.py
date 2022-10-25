from functions.compile_data_tracks_function import *

import ntpath
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind

treatment1 = str(sys.argv[1])

treatment2 = str(sys.argv[2])

min_track_length = int(sys.argv[3])

region1 = str(sys.argv[4])

region2 = str(sys.argv[5])

tracks_region1, tracks_geo_region1, region1_cells, region1_endpointcells = compile_data_tracks(treatment1, min_track_length, region1)

tracks_region2, tracks_geo_region2, region2_cells, region2_endpointcells = compile_data_tracks(treatment2, min_track_length, region2)


#Plot histogram for D/T
plt.hist(region1_endpointcells['DoverT'],bins=30)
plt.title('D/T for {}'.format(region1))
plt.xlabel('D/T')
plt.ylabel('counts')
plt.savefig('histograms/DT_{}.png'.format(region1))

plt.hist(region2_endpointcells['DoverT'],bins=30)
plt.title('D/T for {}'.format(region2))
plt.xlabel('D/T')
plt.ylabel('counts')
plt.savefig('histograms/DT_{}.png'.format(region2))


#Plot histogram for velocity
v_region1 = []
for i in range(len(tracks_geo_region1)):
  v_region1.append(tracks_geo_region1[i]['v'].dropna().tolist())

v_region1 = np.concatenate(v_region1).ravel()

plt.hist(v_region1,bins=50)
plt.title('velocity {}'.format(region1))
plt.savefig('histograms/vel_hist_{}.png'.format(region1))

v_region2 = []
for i in range(len(tracks_geo_region2)):
  v_region2.append(tracks_geo_region2[i]['v'].dropna().tolist())

v_region2 = np.concatenate(v_region2).ravel()

plt.hist(v_region2,bins=50)
plt.title('velocity {}'.format(region2))
plt.savefig('histograms/vel_hist_{}.png'.format(region2))

#Plot histogram for vx and vy
vx_region1 = []
vy_region1 = []
for i in range(len(tracks_geo_region1)):
  vx_region1.append(tracks_geo_region1[i]['vx'].dropna().tolist())
  vy_region1.append(tracks_geo_region1[i]['vy'].dropna().tolist())

vx_region1 = np.concatenate(vx_region1).ravel()
vy_region1 = np.concatenate(vy_region1).ravel()

plt.hist(vx_region1,bins=50)
plt.title('vx {}'.format(region1))
plt.savefig('histograms/velx_hist_{}.png'.format(region1))

plt.hist(vy_region1,bins=50)
plt.title('vy {}'.format(region1))
plt.savefig('histograms/vely_hist_{}.png'.format(region1))

vx_region2 = []
vy_region2 = []
for i in range(len(tracks_geo_region2)):
  vx_region2.append(tracks_geo_region2[i]['vx'].dropna().tolist())
  vy_region2.append(tracks_geo_region2[i]['vy'].dropna().tolist())

vx_region2 = np.concatenate(vx_region2).ravel()
vy_region2 = np.concatenate(vy_region2).ravel()

plt.hist(vx_region2,bins=50)
plt.title('vx {}'.format(region2))
plt.savefig('histograms/velx_hist_{}.png'.format(region2))

plt.hist(vy_region2,bins=50)
plt.title('vy {}'.format(region2))
plt.savefig('histograms/vely_hist_{}.png'.format(region2))

#Plot histogram for dx and dy
dx_region1 = []
dy_region1 = []
for i in range(len(tracks_geo_region1)):
  dx_region1.append(tracks_geo_region1[i]['dx'].dropna())
  dy_region1.append(tracks_geo_region1[i]['dy'].dropna())

dx_region1 = np.concatenate(dx_region1).ravel()
dy_region1 = np.concatenate(dy_region1).ravel()

plt.hist(dx_region1,bins=50)
plt.title('dx {}'.format(region1))
plt.savefig('histograms/dx_hist_{}.png'.format(region1))

plt.hist(dy_region1,bins=50)
plt.title('dy {}'.format(region1))
plt.savefig('histograms/dy_hist_{}.png'.format(region1))

dx_region2 = []
dy_region2 = []
for i in range(len(tracks_geo_region2)):
  dx_region2.append(tracks_geo_region2[i]['dx'].dropna())
  dy_region2.append(tracks_geo_region2[i]['dy'].dropna())

dx_region2 = np.concatenate(dx_region2).ravel()
dy_region2 = np.concatenate(dy_region2).ravel()

plt.hist(dx_region2,bins=50)
plt.title('dx {}'.format(region2))
plt.savefig('histograms/dx_hist_{}.png'.format(region2))

plt.hist(dy_region2,bins=50)
plt.title('dy {}'.format(region2))
plt.savefig('histograms/dy_hist_{}.png'.format(region2))

#make boxplots
data_bp = {'{}'.format(region1):v_region1, '{}'.format(region2):v_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("Velocity")
tstat, pval = ttest_ind(v_region1,v_region2)
plt.text(.2, 20, 'tstatistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))

data_bp = {'{}'.format(region1):vx_region1, '{}'.format(region2):vx_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("Velocity x")
tstat, pval = ttest_ind(vx_region1,vx_region2)
plt.text(.1, 18, 'tstatistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))

data_bp = {'{}'.format(region1):vy_region1, '{}'.format(region2):vy_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("Velocity y")
tstat, pval = ttest_ind(vy_region1,vy_region2)
plt.text(.1, 12, 'tstatistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))

data_bp = {'{}'.format(region1):dx_region1, '{}'.format(region2):dx_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("dx")
tstat, pval = ttest_ind(dx_region1,dx_region2)
plt.text(.1, 12, 'tstatistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))

data_bp = {'{}'.format(region1):dy_region1, '{}'.format(region2):dy_region2}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("dy")
tstat, pval = ttest_ind(dy_region1,dy_region2)
plt.text(.1, 12, 'tstatistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))

data_bp = {'{}'.format(region1):region1_endpointcells['speed'], '{}'.format(region2):region2_endpointcells['speed']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("Speed")
tstat, pval = ttest_ind(region1_endpointcells['speed'],region2_endpointcells['speed'])
plt.text(.1, 15, 'tstatistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))

data_bp = {'{}'.format(region1):region1_endpointcells['DoverT'], '{}'.format(region2):region2_endpointcells['DoverT']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("D/T")
tstat, pval = ttest_ind(region1_endpointcells['DoverT'],region2_endpointcells['DoverT'])
plt.text(.1, 0.8, 'tstatistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))

data_bp = {'{}'.format(region1):region1_endpointcells['FMI'], '{}'.format(region2):region2_endpointcells['FMI']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("FMI")
tstat, pval = ttest_ind(region1_endpointcells['FMI'],region2_endpointcells['FMI'])
plt.text(.1, 0.7, 'tstatistic={}, pvalue={}'.format(round(tstat,3),round(pval,5)), fontsize = 8, bbox = dict(facecolor = 'red', alpha = 0.1))