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

treatment1 = str(sys.argv[1])

treatment2 = str(sys.argv[2])

treatment3 = str(sys.argv[3])

min_track_length = int(sys.argv[4])

region1 = str(sys.argv[5])

region2 = str(sys.argv[6])

region3 = str(sys.argv[7])

save_path = str(sys.argv[8])

pixel_size = 1.54

sampling_t = 5 #min per frame

#Compile dictionaries of data paths indexed by corresponding region name
data_paths={}
if treatment1 != 'False':
  data_paths[region1] = treatment1
if treatment2 != 'False':
  data_paths[region2] = treatment2
if treatment3 != 'False':
  data_paths[region3] = treatment3

#Compile dictionaries for motion metrics calculated from compile_data_tracks function
tracks_region = {}
tracks_geo_region = {}
region_cells = {}
region_endpointcells = {}
for region in data_paths:
  region_data_path = data_paths[region]
  tracks_region[region], tracks_geo_region[region], region_cells[region], region_endpointcells[region] = compile_data_tracks(region_data_path, min_track_length, region, pixel_size)

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

#where acf figures are saved
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
for region in tracks_geo_region:
  absskews[region], solidities[region] = ExtractAbsskewSolidity(tracks_geo_region[region])

#Plot histograms for abs-skew and solidity
for region in absskews:
  HistogramPlot(absskews[region], region, 'absskew', 50, figure_path, file_lines)
  HistogramPlot(solidities[region], region, 'solidity', 50, figure_path, file_lines)


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
region_list = list(tracks_geo_region.keys())

BoxplotPlot(region_list, solidities, 'solidity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(region_list, velocities, 'velocity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(region_list, x_velocities, 'x_velocity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(region_list, y_velocities, 'y_velocity', figure_path, file_lines, pval_file_lines)

BoxplotPlot(region_list, absskews, 'absskew', figure_path, file_lines, pval_file_lines)

BoxplotPlot(region_list, dx, 'dx', figure_path, file_lines, pval_file_lines)

BoxplotPlot(region_list, dy, 'dy', figure_path, file_lines, pval_file_lines)

speeds = {}
for region in region_endpointcells:
  speeds[region] = region_endpointcells[region]['speed']/sampling_t

BoxplotPlot(region_list, speeds, 'speed', figure_path, file_lines, pval_file_lines)

DoverT = {}
for region in region_endpointcells:
  DoverT[region] = region_endpointcells[region]['DoverT']

BoxplotPlot(region_list, DoverT, 'DoverT', figure_path, file_lines, pval_file_lines)

FMI = {}
for region in region_endpointcells:
  FMI[region] = region_endpointcells[region]['FMI']

BoxplotPlot(region_list, FMI, 'FMI', figure_path, file_lines, pval_file_lines)


#track lengths of data
PlotTrackLengthHist(tracks_geo_region, figure_path, file_lines)


#Trajectories
PlotTrajectories(tracks_region, figure_path, file_lines)


#write lines to text file 
hist_boxplot_figs.writelines(file_lines)
hist_boxplot_figs.close() 

pval_tval_txt.writelines(pval_file_lines)
pval_tval_txt.close() 