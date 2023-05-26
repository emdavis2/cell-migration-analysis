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

treatment = str(sys.argv[1])

min_track_length = int(sys.argv[2])

region = str(sys.argv[3])

pixel_size = 1.54

tracks_region, tracks_geo_region, region_cells, region_endpointcells = compile_data_tracks(treatment, min_track_length, region, pixel_size)

if region == 'stiff':
    region_name = 'gel'
else:
    region_name = region

sampling_t = 5 #min per frame

#clears out sentinel file if it exists
open('sentinels/binned_histogram_boxplot_{}.txt'.format(region),'w').close()
#create new sentinel file to write to
hist_boxplot_figs = open('sentinels/binned_histogram_boxplot_{}.txt'.format(region),'w')
file_lines = []

#clears out txt file if it exists
open('figures/binned_histogram_boxplot/pvals_and_tvals_{}.txt'.format(region),'w').close()
#create new txt file to write to
pval_tval_txt = open('figures/binned_histogram_boxplot/pvals_and_tvals_{}.txt'.format(region),'w')
pval_file_lines = []

#bin based on D/T and add column to region_endpointcells dataframe that has label associated with bin
values, bins = np.histogram(region_endpointcells['DoverT'],bins=3)
labels = ['C', 'B', 'A']
region_endpointcells['grade'] = pd.cut(x = region_endpointcells['DoverT'], bins = bins, labels = labels, include_lowest = True)

#filter based on D/T
region_endpointcells_A = region_endpointcells[region_endpointcells['grade'] == 'A'] 
region_endpointcells_B = region_endpointcells[region_endpointcells['grade'] == 'B'] 
region_endpointcells_C = region_endpointcells[region_endpointcells['grade'] == 'C'] 

#names of experiments in each bin
region_experiments_A = []
region_experiments_B = []
region_experiments_C = []
for i in range(len(region_endpointcells_A)):
    key_name = str(region_endpointcells_A.iloc[i]['experiment']) + '/' + str(region_endpointcells_A.iloc[i]['track_id'])
    region_experiments_A.append(key_name)
for i in range(len(region_endpointcells_B)):
    key_name = str(region_endpointcells_B.iloc[i]['experiment']) + '/' + str(region_endpointcells_B.iloc[i]['track_id'])
    region_experiments_B.append(key_name)
for i in range(len(region_endpointcells_C)):
    key_name = str(region_endpointcells_C.iloc[i]['experiment']) + '/' + str(region_endpointcells_C.iloc[i]['track_id'])
    region_experiments_C.append(key_name)

tracks_geo_region_A = []
tracks_region_A = []
tracks_geo_region_B = []
tracks_region_B = []
tracks_geo_region_C = []
tracks_region_C = []
for i in range(len(tracks_geo_region)):
    trial_df = str(tracks_geo_region[i]['experiment'][0]) + '/' + str(tracks_geo_region[i]['track_id'][0])
    if trial_df in region_experiments_A:
        tracks_geo_region_A.append(tracks_geo_region[i])
        tracks_region_A.append(tracks_region[i])
    elif trial_df in region_experiments_B:
        tracks_geo_region_B.append(tracks_geo_region[i])
        tracks_region_B.append(tracks_region[i])
    elif trial_df in region_experiments_C:
        tracks_geo_region_C.append(tracks_geo_region[i])
        tracks_region_C.append(tracks_region[i])
  

#Plot histogram for dx and dy
dx_region_A = []
dy_region_A = []
for i in range(len(tracks_geo_region_A)):
  dx_region_A.append(tracks_geo_region_A[i]['dx'].dropna())
  dy_region_A.append(tracks_geo_region_A[i]['dy'].dropna())

dx_region_A = (np.concatenate(dx_region_A).ravel())/sampling_t
dy_region_A = (np.concatenate(dy_region_A).ravel())/sampling_t

dx_dy_region_A = np.concatenate((dx_region_A,dy_region_A))

plt.hist(dx_dy_region_A,bins=50)
plt.title(r'dx dy ($\mu m$/min) {} D/T min:{} max:{}'.format(region, np.round(bins[2],3), np.round(bins[3],3)))
plt.savefig('figures/binned_histogram_boxplot/dxdy_hist_{}_{}.png'.format(region, 'A'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/dxdy_hist_{}_{}.png \n'.format(region, 'A'))

dx_region_B = []
dy_region_B = []
for i in range(len(tracks_geo_region_B)):
  dx_region_B.append(tracks_geo_region_B[i]['dx'].dropna())
  dy_region_B.append(tracks_geo_region_B[i]['dy'].dropna())

dx_region_B = (np.concatenate(dx_region_B).ravel())/sampling_t
dy_region_B = (np.concatenate(dy_region_B).ravel())/sampling_t

dx_dy_region_B = np.concatenate((dx_region_B,dy_region_B))

plt.hist(dx_dy_region_B,bins=50)
plt.title(r'dx dy ($\mu m$/min) {} D/T min:{} max:{}'.format(region, np.round(bins[1],3), np.round(bins[2],3)))
plt.savefig('figures/binned_histogram_boxplot/dxdy_hist_{}_{}.png'.format(region, 'B'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/dxdy_hist_{}_{}.png \n'.format(region, 'B'))

dx_region_C = []
dy_region_C = []
for i in range(len(tracks_geo_region_C)):
  dx_region_C.append(tracks_geo_region_C[i]['dx'].dropna())
  dy_region_C.append(tracks_geo_region_C[i]['dy'].dropna())

dx_region_C = (np.concatenate(dx_region_C).ravel())/sampling_t
dy_region_C = (np.concatenate(dy_region_C).ravel())/sampling_t

dx_dy_region_C = np.concatenate((dx_region_C,dy_region_C))

plt.hist(dx_dy_region_C,bins=50)
plt.title(r'dx dy ($\mu m$/min) {} D/T min:{} max:{}'.format(region, np.round(bins[0],3), np.round(bins[1],3)))
plt.savefig('figures/binned_histogram_boxplot/dxdy_hist_{}_{}.png'.format(region, 'C'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/dxdy_hist_{}_{}.png \n'.format(region, 'C'))


#Get solidity
solidity_region_A = []
for i in range(len(tracks_geo_region_A)):
  solidity_region_A.append(tracks_geo_region_A[i]['solidity'].dropna().tolist())

solidity_region_A = np.concatenate(solidity_region_A).ravel()

solidity_region_B = []
for i in range(len(tracks_geo_region_B)):
  solidity_region_B.append(tracks_geo_region_B[i]['solidity'].dropna().tolist())

solidity_region_B = np.concatenate(solidity_region_B).ravel()

solidity_region_C = []
for i in range(len(tracks_geo_region_C)):
  solidity_region_C.append(tracks_geo_region_C[i]['solidity'].dropna().tolist())

solidity_region_C = np.concatenate(solidity_region_C).ravel()


#make boxplots
data_bp = {'D/T {} - {}'.format(np.round(bins[0],3), np.round(bins[1],3)):solidity_region_C, 'D/T {} - {}'.format(np.round(bins[1],3), np.round(bins[2],3)):solidity_region_B, 'D/T {} - {}'.format(np.round(bins[2],3), np.round(bins[3],3)):solidity_region_A}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("Solidity on {}".format(region))
tstatAB, pvalAB = f_oneway(solidity_region_A,solidity_region_B)
tstatAC, pvalAC = f_oneway(solidity_region_A,solidity_region_C)
tstatBC, pvalBC = f_oneway(solidity_region_B,solidity_region_C)
pval_file_lines.append('\n Solidity \n')
pval_file_lines.append('{} and {}: {} \n'.format('A', 'B', pvalAB))
pval_file_lines.append('{} and {}: {} \n'.format('A', 'C', pvalAC))
pval_file_lines.append('{} and {}: {} \n'.format('B', 'C', pvalBC))
# plt.text(0.1, 0.8, 'pvalue_12={}'.format(pval12), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(0.1, 0.7, 'pvalue_13={}'.format(pval13), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(0.1, 0.6, 'pvalue_23={}'.format(pval23), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/binned_histogram_boxplot/solidity_boxplot_{}.png'.format(region))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/solidity_boxplot_{}.png \n'.format(region))

data_bp = {'D/T {} - {}'.format(np.round(bins[0],3), np.round(bins[1],3)):region_endpointcells_C['speed']/sampling_t, 'D/T {} - {}'.format(np.round(bins[1],3), np.round(bins[2],3)):region_endpointcells_B['speed']/sampling_t, 'D/T {} - {}'.format(np.round(bins[2],3), np.round(bins[3],3)):region_endpointcells_A['speed']/sampling_t}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel(r"Speed ($\mu m$/min) on {}".format(region))
tstatAB, pvalAB = f_oneway(region_endpointcells_A['speed'],region_endpointcells_B['speed'])
tstatAC, pvalAC = f_oneway(region_endpointcells_A['speed'],region_endpointcells_C['speed'])
tstatBC, pvalBC = f_oneway(region_endpointcells_B['speed'],region_endpointcells_C['speed'])
pval_file_lines.append('\n Speed \n')
pval_file_lines.append('{} and {}: {} \n'.format('A', 'B', pvalAB))
pval_file_lines.append('{} and {}: {} \n'.format('A', 'C', pvalAC))
pval_file_lines.append('{} and {}: {} \n'.format('B', 'C', pvalBC))
# plt.text(.1, 1.5, 'pvalue_12={}'.format(pval12), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(.1, 1.2, 'pvalue_13={}'.format(pval13), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(.1, 1, 'pvalue_23={}'.format(pval23), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/binned_histogram_boxplot/speed_boxplot_{}.png'.format(region))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/speed_boxplot_{}.png \n'.format(region))

data_bp = {'D/T {} - {}'.format(np.round(bins[0],3), np.round(bins[1],3)):region_endpointcells_C['DoverT'], 'D/T {} - {}'.format(np.round(bins[1],3), np.round(bins[2],3)):region_endpointcells_B['DoverT'], 'D/T {} - {}'.format(np.round(bins[2],3), np.round(bins[3],3)):region_endpointcells_A['DoverT']}
data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
sns.boxplot(data=data_boxplot)
plt.xlabel("Source")
plt.ylabel("D/T on {}".format(region))
tstatAB, pvalAB = f_oneway(region_endpointcells_A['DoverT'],region_endpointcells_B['DoverT'])
tstatAC, pvalAC = f_oneway(region_endpointcells_A['DoverT'],region_endpointcells_C['DoverT'])
tstatBC, pvalBC = f_oneway(region_endpointcells_B['DoverT'],region_endpointcells_C['DoverT'])
pval_file_lines.append('\n D/T \n')
pval_file_lines.append('{} and {}: {} \n'.format('A', 'B', pvalAB))
pval_file_lines.append('{} and {}: {} \n'.format('A', 'C', pvalAC))
pval_file_lines.append('{} and {}: {} \n'.format('B', 'C', pvalBC))
# plt.text(.1, 0.4, 'pvalue_12={}'.format(pval12), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(.1, 0.3, 'pvalue_13={}'.format(pval13), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
# plt.text(.1, 0.2, 'pvalue_23={}'.format(pval23), fontsize = 12, bbox = dict(facecolor = 'red', alpha = 0.1))
plt.savefig('figures/binned_histogram_boxplot/DoverT_boxplot_{}.png'.format(region))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/DoverT_boxplot_{}.png \n'.format(region))


#track lengths of data
lengths = [len(tracks_geo_region_A[i]) for i in range(len(tracks_geo_region_A)) ]
num_tracks = len(tracks_geo_region_A)
plt.hist(lengths)
plt.xlabel('track length (number of frames)')
plt.ylabel('counts')
plt.title('Track lengths for {} data D/T min:{} max:{} ({} tracks)'.format(region_name, np.round(bins[2],3), np.round(bins[3],3), num_tracks))
plt.savefig('figures/binned_histogram_boxplot/track_lengths_data_{}_{}.png'.format(region, 'A'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/track_lengths_data_{}_{}.png \n'.format(region, 'A'))

#track lengths of data
lengths = [len(tracks_geo_region_B[i]) for i in range(len(tracks_geo_region_B)) ]
num_tracks = len(tracks_geo_region_B)
plt.hist(lengths)
plt.xlabel('track length (number of frames)')
plt.ylabel('counts')
plt.title('Track lengths for {} data D/T min:{} max:{} ({} tracks)'.format(region_name, np.round(bins[1],3), np.round(bins[2],3), num_tracks))
plt.savefig('figures/binned_histogram_boxplot/track_lengths_data_{}_{}.png'.format(region, 'B'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/track_lengths_data_{}_{}.png \n'.format(region, 'B'))

#track lengths of data
lengths = [len(tracks_geo_region_C[i]) for i in range(len(tracks_geo_region_C)) ]
num_tracks = len(tracks_geo_region_C)
plt.hist(lengths)
plt.xlabel('track length (number of frames)')
plt.ylabel('counts')
plt.title('Track lengths for {} data D/T min:{} max:{} ({} tracks)'.format(region_name, np.round(bins[0],3), np.round(bins[1],3), num_tracks))
plt.savefig('figures/binned_histogram_boxplot/track_lengths_data_{}_{}.png'.format(region, 'C'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/track_lengths_data_{}_{}.png \n'.format(region, 'C'))

#Trajectories
for df in tracks_region_A:
    plt.plot(df[:,0],df[:,1])
    #x_data[i] = df[:,0]
    #y_data[i] = df[:,1]
    #i += 1
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for {} Data D/T min:{} max:{}'.format(region_name, np.round(bins[2],3), np.round(bins[3],3)))
plt.savefig('figures/binned_histogram_boxplot/trajectories_data_{}_{}.png'.format(region, 'A'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/trajectories_data_{}.png \n'.format(region, 'A'))

#Trajectories
for df in tracks_region_B:
    plt.plot(df[:,0],df[:,1])
    #x_data[i] = df[:,0]
    #y_data[i] = df[:,1]
    #i += 1
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for {} Data D/T min:{} max:{}'.format(region_name, np.round(bins[1],3), np.round(bins[2],3)))
plt.savefig('figures/binned_histogram_boxplot/trajectories_data_{}_{}.png'.format(region, 'B'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/trajectories_data_{}.png \n'.format(region, 'B'))

#Trajectories
for df in tracks_region_C:
    plt.plot(df[:,0],df[:,1])
    #x_data[i] = df[:,0]
    #y_data[i] = df[:,1]
    #i += 1
plt.xlabel(('x position ($\mu$m)'))
plt.ylabel(('y position ($\mu$m)'))
plt.title('Trajectories for {} Data D/T min:{} max:{}'.format(region_name, np.round(bins[0],3), np.round(bins[1],3)))
plt.savefig('figures/binned_histogram_boxplot/trajectories_data_{}_{}.png'.format(region, 'C'))
plt.clf()
file_lines.append('figures/binned_histogram_boxplot/trajectories_data_{}.png \n'.format(region, 'C'))

#write lines to text file 
hist_boxplot_figs.writelines(file_lines)
hist_boxplot_figs.close() 

#write lines to text file 
pval_tval_txt.writelines(pval_file_lines)
pval_tval_txt.close() 