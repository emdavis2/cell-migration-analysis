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

pixel_size = 1.54 #pixels/micron

tracks_region1, tracks_geo_region1, region1_cells, region1_endpointcells = compile_data_tracks(treatment1, min_track_length, region1, pixel_size)

tracks_region2, tracks_geo_region2, region2_cells, region2_endpointcells = compile_data_tracks(treatment2, min_track_length, region2, pixel_size)

tracks_region3, tracks_geo_region3, region3_cells, region3_endpointcells = compile_data_tracks(treatment3, min_track_length, region3, pixel_size)


#clears out sentinel file if it exists
open('sentinels/cellshape_histogram.txt','w').close()
#create new sentinel file to write to
hist_boxplot_figs = open('sentinels/cellshape_histogram.txt','w')
file_lines = [] 

var_area_region1 = []
area_region1 = []
for df in tracks_geo_region1:
  area = df['area'].iloc[:,0]/(pixel_size**2)
  area_region1.append(area[:min_track_length])
  var_area_region1.append(np.var(area))
avg_area_region1 = np.average(area_region1, axis=0)
plt.title('Area over time (over track length)')
plt.xlabel('time (5 minutes)')
plt.ylabel('area (microns squared)')
plt.plot(area_region1)
plt.savefig('figures/cellshape_histogram/area_overtime_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/cellshape_histogram/area_overtime_{}.png \n'.format(region1))
plt.title('Average area over time (over track length of all cells)')
plt.xlabel('time (5 minutes)')
plt.ylabel('area (microns squared)')
plt.plot(avg_area_region1)
plt.savefig('figures/cellshape_histogram/avg_area_overtime_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/cellshape_histogram/avg_area_overtime_{}.png \n'.format(region1))
# plt.title('Variance of cell area over track in {}'.format(region1))
# plt.hist(var_area_region1, bins=30)
# plt.savefig('figures/cellshape_histogram/var_area_hist_{}.png'.format(region1))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/var_area_hist_{}.png \n'.format(region1))
# plt.title('Cell area over track in {}'.format(region1))
# plt.hist(np.concatenate(area_region1).ravel(), bins=30)
# plt.savefig('figures/cellshape_histogram/area_hist_{}.png'.format(region1))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/area_hist_{}.png \n'.format(region1))

var_area_region2 = []
area_region2 = []
for df in tracks_geo_region2:
  area = df['area'].iloc[:,0]/(pixel_size**2)
  area_region2.append(area[:min_track_length])
  var_area_region2.append(np.var(area))
avg_area_region2 = np.average(area_region2, axis=0)
plt.title('Area over time (over track length)')
plt.xlabel('time (5 minutes)')
plt.ylabel('area (microns squared)')
plt.plot(area_region2)
plt.savefig('figures/cellshape_histogram/area_overtime_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/cellshape_histogram/area_overtime_{}.png \n'.format(region2))
plt.title('Average area over time (over track length of all cells)')
plt.xlabel('time (5 minutes)')
plt.ylabel('area (microns squared)')
plt.plot(avg_area_region2)
plt.savefig('figures/cellshape_histogram/avg_area_overtime_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/cellshape_histogram/avg_area_overtime_{}.png \n'.format(region2))
# plt.title('Variance of cell area over track on {}'.format(region2))
# plt.hist(var_area_region2, bins=30)
# plt.savefig('figures/cellshape_histogram/var_area_hist_{}.png'.format(region2))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/var_area_hist_{}.png \n'.format(region2))
# plt.title('Cell area over track in {}'.format(region2))
# plt.hist(np.concatenate(area_region2).ravel(), bins=30)
# plt.savefig('figures/cellshape_histogram/area_hist_{}.png'.format(region2))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/area_hist_{}.png \n'.format(region2))

var_area_region3 = []
area_region3 = []
for df in tracks_geo_region3:
  area = df['area'].iloc[:,0]/(pixel_size**2)
  area_region3.append(area[:min_track_length])
  var_area_region3.append(np.var(area))
avg_area_region3 = np.average(area_region3, axis=0)
plt.title('Area over time (over track length)')
plt.xlabel('time (5 minutes)')
plt.ylabel('area (microns squared)')
plt.plot(area_region3)
plt.savefig('figures/cellshape_histogram/area_overtime_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/cellshape_histogram/area_overtime_{}.png \n'.format(region3))
plt.title('Average area over time (over track length of all cells)')
plt.xlabel('time (5 minutes)')
plt.ylabel('area (microns squared)')
plt.plot(avg_area_region3)
plt.savefig('figures/cellshape_histogram/avg_area_overtime_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/cellshape_histogram/avg_area_overtime_{}.png \n'.format(region3))
# plt.title('Variance of cell area over track on {}'.format(region3))
# plt.hist(var_area_region3, bins=30)
# plt.savefig('figures/cellshape_histogram/var_area_hist_{}.png'.format(region3))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/var_area_hist_{}.png \n'.format(region3))
# plt.title('Cell area over track in {}'.format(region3))
# plt.hist(np.concatenate(area_region3).ravel(), bins=30)
# plt.savefig('figures/cellshape_histogram/area_hist_{}.png'.format(region3))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/area_hist_{}.png \n'.format(region3))

var_perim_region1 = []
perim_region1 = []
for df in tracks_geo_region1:
  perim = df['perimeter']/pixel_size
  perim_region1.append(perim[:min_track_length])
  var_perim_region1.append(np.var(perim))
avg_perim_region1 = np.average(perim_region1, axis=0)
plt.title('Perimeter over time (over track length)')
plt.xlabel('time (5 minutes)')
plt.ylabel('perimeter (microns)')
plt.plot(perim_region1)
plt.savefig('figures/cellshape_histogram/perim_overtime_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/cellshape_histogram/perim_overtime_{}.png \n'.format(region1))
plt.title('Average perimeter over time (over track length of all cells)')
plt.xlabel('time (5 minutes)')
plt.ylabel('perimeter (microns)')
plt.plot(avg_perim_region1)
plt.savefig('figures/cellshape_histogram/avg_perim_overtime_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/cellshape_histogram/avg_perim_overtime_{}.png \n'.format(region1))
# plt.title('Variance of cell perimeter over track on {}'.format(region1))
# plt.hist(var_perim_region1, bins=30)
# plt.savefig('figures/cellshape_histogram/var_perim_hist_{}.png'.format(region1))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/var_perim_hist_{}.png \n'.format(region1))
# plt.title('Cell perimeter over track on {}'.format(region1))
# plt.hist(np.concatenate(perim_region1).ravel(), bins=30)
# plt.savefig('figures/cellshape_histogram/perim_hist_{}.png'.format(region1))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/perim_hist_{}.png \n'.format(region1))

var_perim_region2 = []
perim_region2 = []
for df in tracks_geo_region2:
  perim = df['perimeter']/pixel_size
  perim_region2.append(perim[:min_track_length])
  var_perim_region2.append(np.var(perim))
avg_perim_region2 = np.average(perim_region2, axis=0)
plt.title('Perimeter over time (over track length)')
plt.xlabel('time (5 minutes)')
plt.ylabel('perimeter (microns)')
plt.plot(perim_region2)
plt.savefig('figures/cellshape_histogram/perim_overtime_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/cellshape_histogram/perim_overtime_{}.png \n'.format(region2))
plt.title('Average perimeter over time (over track length of all cells)')
plt.xlabel('time (5 minutes)')
plt.ylabel('perimeter (microns)')
plt.plot(avg_perim_region2)
plt.savefig('figures/cellshape_histogram/avg_perim_overtime_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/cellshape_histogram/avg_perim_overtime_{}.png \n'.format(region2))
# plt.title('Variance of cell perimeter over track on {}'.format(region2))
# plt.hist(var_perim_region2, bins=30)
# plt.savefig('figures/cellshape_histogram/perim_hist_{}.png'.format(region2))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/perim_hist_{}.png \n'.format(region2))
# plt.title('Cell perimeter over track on {}'.format(region2))
# plt.hist(np.concatenate(perim_region2).ravel(), bins=30)
# plt.savefig('figures/cellshape_histogram/perim_hist_{}.png'.format(region2))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/perim_hist_{}.png \n'.format(region2))

var_perim_region3 = []
perim_region3 = []
for df in tracks_geo_region3:
  perim = df['perimeter']/pixel_size
  perim_region3.append(perim[:min_track_length])
  var_perim_region3.append(np.var(perim))
avg_perim_region3 = np.average(perim_region3, axis=0)
plt.title('Perimeter over time (over track length)')
plt.xlabel('time (5 minutes)')
plt.ylabel('perimeter (microns)')
plt.plot(perim_region3)
plt.savefig('figures/cellshape_histogram/perim_overtime_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/cellshape_histogram/perim_overtime_{}.png \n'.format(region3))
plt.title('Average perimeter over time (over track length of all cells)')
plt.xlabel('time (5 minutes)')
plt.ylabel('perimeter (microns)')
plt.plot(avg_perim_region3)
plt.savefig('figures/cellshape_histogram/avg_perim_overtime_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/cellshape_histogram/avg_perim_overtime_{}.png \n'.format(region3))
# plt.title('Variance of cell perimeter over track on {}'.format(region3))
# plt.hist(var_perim_region3, bins=30)
# plt.savefig('figures/cellshape_histogram/perim_hist_{}.png'.format(region3))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/perim_hist_{}.png \n'.format(region3))
# plt.title('Cell perimeter over track on {}'.format(region3))
# plt.hist(np.concatenate(perim_region3).ravel(), bins=30)
# plt.savefig('figures/cellshape_histogram/perim_hist_{}.png'.format(region3))
# plt.clf()
# file_lines.append('figures/cellshape_histogram/perim_hist_{}.png \n'.format(region3))

###########PERCENT CHANGE#######################


perchange_area_region1 = []
perchange_perim_region1 = []
for df in tracks_geo_region1:
  area = df['area'].iloc[:,0]/(pixel_size**2)
  perim = df['perimeter']/pixel_size
  max_area = np.max(area)
  min_area = np.min(area)
  perchange_area_region1.append((max_area-min_area)/max_area)
  max_perim = np.max(perim)
  min_perim = np.min(perim)
  perchange_perim_region1.append((max_perim-min_perim)/max_perim)
plt.title('Percent change of cell area over track in {}'.format(region1))
plt.hist(perchange_area_region1)
plt.savefig('figures/cellshape_histogram/perchange_area_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/cellshape_histogram/perchange_area_hist_{}.png \n'.format(region1))
plt.title('Percent change of cell perimeter over track in {}'.format(region1))
plt.hist(perchange_perim_region1)
plt.savefig('figures/cellshape_histogram/perchange_perim_hist_{}.png'.format(region1))
plt.clf()
file_lines.append('figures/cellshape_histogram/perchange_perim_hist_{}.png \n'.format(region1))

perchange_area_region2 = []
perchange_perim_region2 = []
for df in tracks_geo_region2:
  area = df['area'].iloc[:,0]/(pixel_size**2)
  perim = df['perimeter']/pixel_size
  max_area = np.max(area)
  min_area = np.min(area)
  perchange_area_region2.append((max_area-min_area)/max_area)
  max_perim = np.max(perim)
  min_perim = np.min(perim)
  perchange_perim_region2.append((max_perim-min_perim)/max_perim)
plt.title('Percent change of cell area over track in {}'.format(region2))
plt.hist(perchange_area_region2)
plt.savefig('figures/cellshape_histogram/perchange_area_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/cellshape_histogram/perchange_area_hist_{}.png \n'.format(region2))
plt.title('Percent change of cell perimeter over track in {}'.format(region2))
plt.hist(perchange_perim_region2)
plt.savefig('figures/cellshape_histogram/perchange_perim_hist_{}.png'.format(region2))
plt.clf()
file_lines.append('figures/cellshape_histogram/perchange_perim_hist_{}.png \n'.format(region2))

perchange_area_region3 = []
perchange_perim_region3 = []
for df in tracks_geo_region1:
  area = df['area'].iloc[:,0]/(pixel_size**2)
  perim = df['perimeter']/pixel_size
  max_area = np.max(area)
  min_area = np.min(area)
  perchange_area_region3.append((max_area-min_area)/max_area)
  max_perim = np.max(perim)
  min_perim = np.min(perim)
  perchange_perim_region3.append((max_perim-min_perim)/max_perim)
plt.title('Percent change of cell area over track in {}'.format(region3))
plt.hist(perchange_area_region3)
plt.savefig('figures/cellshape_histogram/perchange_area_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/cellshape_histogram/perchange_area_hist_{}.png \n'.format(region3))
plt.title('Percent change of cell perimeter over track in {}'.format(region3))
plt.hist(perchange_perim_region3)
plt.savefig('figures/cellshape_histogram/perchange_perim_hist_{}.png'.format(region3))
plt.clf()
file_lines.append('figures/cellshape_histogram/perchange_perim_hist_{}.png \n'.format(region3))