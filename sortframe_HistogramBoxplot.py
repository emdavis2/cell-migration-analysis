from functions.compile_data_tracks_function import *
from functions.plot_functions import *

import ntpath
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
#from scipy.stats import ttest_ind
from scipy.stats import f_oneway

### Treatment input is the list of paths to the data-the data should be structured such that the pkl files are organized as:
### gel_region > treatment_type/day > pkl files

treatment_input = sys.argv[1]

min_track_length = int(sys.argv[2])

region_input = sys.argv[3]

save_path = str(sys.argv[4])

pixel_size = 1.54

sampling_t = 5 #min per frame

#turn treatment and region into lists of strings
treatment_list = list(map(str, treatment_input.strip('[]').split(',')))
region_list = list(map(str, region_input.strip('[]').split(',')))

#Get list of unique defining names for each data path
treatment_names =[]
for ind,path in enumerate(treatment_list):
  treatment_names.append(region_list[ind]+'_'+ntpath.basename(path))

#Compile dictionaries of data paths indexed by corresponding region name
data_paths={}
for ind,treatment in enumerate(treatment_list):
  if treatment != 'False':
    data_paths[treatment_names[ind]] = treatment

#Compile dictionaries for motion metrics calculated from compile_data_tracks function
tracks_region = {}
tracks_geo_region = {}
region_cells = {}
region_endpointcells = {}
for ind,treat_name in enumerate(data_paths):
  region_data_path = data_paths[treat_name]
  tracks_region[treat_name], tracks_geo_region[treat_name], region_cells[treat_name], region_endpointcells[treat_name] = compile_data_tracks(region_data_path, min_track_length, region_list[ind], pixel_size)

#Determine which tracks are in first half of movie and which are in second half of movie
filtered_tracks_geo_region = {}
filtered_tracks_region = {}
filtered_region_endpointcells = {}
# experiment_names = {}
for treat_name in tracks_geo_region:
  filtered_tracks_geo_region['firsthalf_{}'.format(treat_name)] = []
  filtered_tracks_geo_region['secondhalf_{}'.format(treat_name)] = []
  filtered_tracks_region['firsthalf_{}'.format(treat_name)] = []
  filtered_tracks_region['secondhalf_{}'.format(treat_name)] = []
  filtered_region_endpointcells['firsthalf_{}'.format(treat_name)] = []
  filtered_region_endpointcells['secondhalf_{}'.format(treat_name)] = []
#   experiment_names['firsthalf_{}'.format(treat_name)] = []
#   experiment_names['secondhalf_{}'.format(treat_name)] = []
  for i,df in enumerate(tracks_geo_region[treat_name]):
    if df['frame'][0] < 108:
      filtered_tracks_geo_region['firsthalf_{}'.format(treat_name)].append(df)
      filtered_tracks_region['firsthalf_{}'.format(treat_name)].append(tracks_region[treat_name][i])
      filtered_region_endpointcells['firsthalf_{}'.format(treat_name)] = region_endpointcells[treat_name][region_endpointcells[treat_name]['experiment'] == df['experiment'][0]]
    #   experiment_names['firsthalf_{}'.format(treat_name)].append(df['experiment'][0])
    else:
      filtered_tracks_geo_region['secondhalf_{}'.format(treat_name)].append(df)
      filtered_tracks_region['secondhalf_{}'.format(treat_name)].append(tracks_region[treat_name][i])
      filtered_region_endpointcells['secondhalf_{}'.format(treat_name)] = region_endpointcells[treat_name][region_endpointcells[treat_name]['experiment'] == df['experiment'][0]]
    #   experiment_names['secondhalf_{}'.format(treat_name)].append(df['experiment'][0])



tracks_geo_region = filtered_tracks_geo_region
tracks_region = filtered_tracks_region
region_endpointcells = filtered_region_endpointcells

treatment_names = list(tracks_geo_region.keys())

#check to see if the path exists, if not make the directory
if not os.path.exists('sentinels'):
  os.mkdir('sentinels')

#clears out sentinel file if it exists
open('sentinels/histogram_boxplot.txt','w').close()
#create new sentinel file to write to
hist_boxplot_figs = open('sentinels/histogram_boxplot.txt','w')
file_lines = []

#check to see if the path exists, if not make the directory
if not os.path.exists(save_path):
  os.mkdir(save_path)

#where figures are saved
figure_path = save_path + '/histogram_boxplot'

#check to see if the path exists, if not make the directory
if not os.path.exists(figure_path):
  os.mkdir(figure_path)

#clears out txt file if it exists
open('{}/pvals_and_tvals.txt'.format(figure_path),'w').close()
#create new txt file to write to
pval_tval_txt = open('{}/pvals_and_tvals.txt'.format(figure_path),'w')
pval_file_lines = []

#Plot histogram for D/T
for region in region_endpointcells:
  HistogramPlot(region_endpointcells[region]['DoverT'], region, 'DoverT', 30, figure_path, file_lines)


#Plot histogram for speed
for region in region_endpointcells:
  HistogramPlot(region_endpointcells[region]['speed'], region, 'speed', 30, figure_path, file_lines)


#Plot histogram for FMI
for region in region_endpointcells:
  HistogramPlot(region_endpointcells[region]['FMI'], region, 'FMI', 30, figure_path, file_lines)

#Compile a dictionary of velocities indexed by region
velocities = {}
x_velocities = {}
y_velocities = {}
for region in tracks_geo_region:
  velocities[region], x_velocities[region], y_velocities[region] = ExtractVelocity(tracks_geo_region[region], sampling_t)

#Plot histogram for velocity
for region in velocities:
  HistogramPlot(velocities[region], region, 'velocity', 50, figure_path, file_lines)


#Plot histogram for vx and vy
for region in x_velocities:
  HistogramPlot(x_velocities[region], region, 'x_velocity', 50, figure_path, file_lines)
  HistogramPlot(y_velocities[region], region, 'y_velocity', 50, figure_path, file_lines)


#Compile a dictionary of abs-skew and solidity indexed by region
absskews = {}
solidities = {}
eccentricities = {}
areas = {}
for region in tracks_geo_region:
  absskews[region], solidities[region], eccentricities[region], areas[region] = ExtractAbsskewSolidityEccenArea(tracks_geo_region[region], pixel_size)

#Plot histograms for abs-skew and solidity
for region in absskews:
  HistogramPlot(absskews[region], region, 'absskew', 50, figure_path, file_lines)
  HistogramPlot(solidities[region], region, 'solidity', 50, figure_path, file_lines)
  HistogramPlot(eccentricities[region], region, 'eccentricity', 50, figure_path, file_lines)


#Compile a dictionary of dx and dy indexed by regions
dx = {}
dy = {}
dx_dy = {}
for region in tracks_geo_region:
  dx_region = []
  dy_region = []
  for i in range(len(tracks_geo_region[region])):
    dx_region.append(tracks_geo_region[region][i]['dx'].dropna())
    dy_region.append(tracks_geo_region[region][i]['dy'].dropna())
  dx[region] = (np.concatenate(dx_region).ravel())/sampling_t
  dy[region] = (np.concatenate(dy_region).ravel())/sampling_t
  dx_dy[region] = np.concatenate((dx[region], dy[region]))

# Plot histogram for dx and dy
for region in dx:
  HistogramPlot(dx[region], region, 'dx', 50, figure_path, file_lines)
  HistogramPlot(dy[region], region, 'dy', 50, figure_path, file_lines)
  HistogramPlot(dx_dy[region], region, 'dxdy', 50, figure_path, file_lines)

#make boxplots

BoxplotPlot(treatment_names, solidities, 'solidity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(treatment_names, eccentricities, 'eccentricity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(treatment_names, areas, 'area', figure_path, file_lines, pval_file_lines)

BoxplotPlot(treatment_names, velocities, 'velocity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(treatment_names, x_velocities, 'x_velocity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(treatment_names, y_velocities, 'y_velocity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(treatment_names, absskews, 'absskew', figure_path, file_lines, pval_file_lines)

BoxplotPlot(treatment_names, dx, 'dx', figure_path, file_lines, pval_file_lines)

BoxplotPlot(treatment_names, dy, 'dy', figure_path, file_lines, pval_file_lines)

speeds = {}
for region in region_endpointcells:
  speeds[region] = region_endpointcells[region]['speed']/sampling_t

BoxplotPlot(treatment_names, speeds, 'speed', figure_path, file_lines, pval_file_lines)

DoverT = {}
for region in region_endpointcells:
  DoverT[region] = region_endpointcells[region]['DoverT']

BoxplotPlot(treatment_names, DoverT, 'DoverT', figure_path, file_lines, pval_file_lines)

FMI = {}
for region in region_endpointcells:
  FMI[region] = region_endpointcells[region]['FMI']

BoxplotPlot(treatment_names, FMI, 'FMI', figure_path, file_lines, pval_file_lines)


#track lengths of data
PlotTrackLengthHist(tracks_geo_region, figure_path, file_lines)


#Trajectories
PlotTrajectories(tracks_region, figure_path, file_lines)


#write lines to text file 
hist_boxplot_figs.writelines(file_lines)
hist_boxplot_figs.close() 

pval_tval_txt.writelines(pval_file_lines)
pval_tval_txt.close() 