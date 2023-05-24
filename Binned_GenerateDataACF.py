from functions.acf_functions import *
from functions.compile_data_tracks_function import *

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

treatment = str(sys.argv[1])

min_track_length = int(sys.argv[2])

region = str(sys.argv[3])

pixel_size = 1.54

min_persistence = 0.6

tracks_region, tracks_geo_region, region_cells, region_endpointcells = compile_data_tracks(treatment, min_track_length, region, pixel_size)

if region == 'stiff':
  region_name = 'gel'
else:
  region_name = region

#clears out sentinel file if it exists
open('sentinels/binned_ACF_figures_{}.txt'.format(region),'w').close()
#create new sentinel file to write to
acf_fig_region = open('sentinels/binned_ACF_figures_{}.txt'.format(region),'w')
file_lines = []

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


#autocorrelation velocity 
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region_A:
  track=df
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
  #combined = combined.iloc[1::2]
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_A) #Nposlagtotal 

std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel('lag (5 min)')
plt.title("Autocorrelaton velocity {} for D/T {} - {}".format(region_name, np.round(bins[2],3), np.round(bins[3],3)))
plt.savefig('figures/binned_acf_figures/{}_velocity_acf_avg_{}.png'.format(region, 'A'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_velocity_acf_avg_{}.png \n'.format(region, 'A'))

#autocorrelation velocity 
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region_B:
  track=df
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
  #combined = combined.iloc[1::2]
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_B) #Nposlagtotal 

std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel('lag (5 min)')
plt.title("Autocorrelaton velocity {} for D/T {} - {}".format(region_name, np.round(bins[1],3), np.round(bins[2],3)))
plt.savefig('figures/binned_acf_figures/{}_velocity_acf_avg_{}.png'.format(region, 'B'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_velocity_acf_avg_{}.png \n'.format(region, 'B'))

#autocorrelation velocity 
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region_C:
  track=df
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
  #combined = combined.iloc[1::2]
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_C) #Nposlagtotal 

std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel('lag (5 min)')
plt.title("Autocorrelaton velocity {} for D/T {} - {}".format(region_name, np.round(bins[2],3), np.round(bins[3],3)))
plt.savefig('figures/binned_acf_figures/{}_velocity_acf_avg_{}.png'.format(region, 'C'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_velocity_acf_avg_{}.png \n'.format(region, 'C'))

#cross correlation velocity vector and polarity vector
poslagaverage = np.zeros(300)
neglagaverage = np.zeros(300)
#Nposlagtotal = np.zeros(300)
all_ac_pos = []
all_ac_neg = []
for df in tracks_geo_region_A:
  track=df
  cospol = list( track[['abs-skew']].iloc[:,0] * np.cos(track[['polarity_angle']].iloc[:,0]) )
  sinpol = list( track[['abs-skew']].iloc[:,0] * np.sin(track[['polarity_angle']].iloc[:,0]) ) 
  normvect = pd.DataFrame( {'cospol': cospol , 'sinpol':sinpol })
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  neglagsmean[np.isnan(neglagsmean)] = 0
  all_ac_pos.append(poslagsmean)
  all_ac_neg.append(neglagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
  neglagaverage[0:len(neglagsmean)] += neglagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_A) #Nposlagtotal 
neglagaverage /= len(tracks_geo_region_A)

std_err_pos = np.std(all_ac_pos,axis=0,ddof=1)/np.sqrt(np.shape(all_ac_pos)[0])
std_err_neg = np.std(all_ac_neg,axis=0,ddof=1)/np.sqrt(np.shape(all_ac_neg)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
#plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err_pos,label='pos lag')
plt.errorbar(np.arange(0,min_track_length-4),neglagaverage[0:min_track_length-4],yerr=std_err_neg,label='neg lag')
plt.xlabel('lag (10 min)')
plt.title("Cross correlation velocity and polarity vectors {} for D/T {} - {}".format(region_name, np.round(bins[2],3), np.round(bins[3],3)))
plt.savefig('figures/binned_acf_figures/{}_vel_pol_crosscorr_avg_{}.png'.format(region, 'A'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_vel_pol_crosscorr_avg_{}.png \n'.format(region, 'A'))

#cross correlation velocity vector and polarity vector
poslagaverage = np.zeros(300)
neglagaverage = np.zeros(300)
#Nposlagtotal = np.zeros(300)
all_ac_pos = []
all_ac_neg = []
for df in tracks_geo_region_B:
  track=df
  cospol = list( track[['abs-skew']].iloc[:,0] * np.cos(track[['polarity_angle']].iloc[:,0]) )
  sinpol = list( track[['abs-skew']].iloc[:,0] * np.sin(track[['polarity_angle']].iloc[:,0]) ) 
  normvect = pd.DataFrame( {'cospol': cospol , 'sinpol':sinpol })
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  neglagsmean[np.isnan(neglagsmean)] = 0
  all_ac_pos.append(poslagsmean)
  all_ac_neg.append(neglagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
  neglagaverage[0:len(neglagsmean)] += neglagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_B) #Nposlagtotal 
neglagaverage /= len(tracks_geo_region_B)

std_err_pos = np.std(all_ac_pos,axis=0,ddof=1)/np.sqrt(np.shape(all_ac_pos)[0])
std_err_neg = np.std(all_ac_neg,axis=0,ddof=1)/np.sqrt(np.shape(all_ac_neg)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
#plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err_pos,label='pos lag')
plt.errorbar(np.arange(0,min_track_length-4),neglagaverage[0:min_track_length-4],yerr=std_err_neg,label='neg lag')
plt.xlabel('lag (10 min)')
plt.title("Cross correlation velocity and polarity vectors {} for D/T {} - {}".format(region_name, np.round(bins[1],3), np.round(bins[2],3)))
plt.savefig('figures/binned_acf_figures/{}_vel_pol_crosscorr_avg_{}.png'.format(region, 'B'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_vel_pol_crosscorr_avg_{}.png \n'.format(region, 'B'))

#cross correlation velocity vector and polarity vector
poslagaverage = np.zeros(300)
neglagaverage = np.zeros(300)
#Nposlagtotal = np.zeros(300)
all_ac_pos = []
all_ac_neg = []
for df in tracks_geo_region_C:
  track=df
  cospol = list( track[['abs-skew']].iloc[:,0] * np.cos(track[['polarity_angle']].iloc[:,0]) )
  sinpol = list( track[['abs-skew']].iloc[:,0] * np.sin(track[['polarity_angle']].iloc[:,0]) ) 
  normvect = pd.DataFrame( {'cospol': cospol , 'sinpol':sinpol })
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  neglagsmean[np.isnan(neglagsmean)] = 0
  all_ac_pos.append(poslagsmean)
  all_ac_neg.append(neglagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
  neglagaverage[0:len(neglagsmean)] += neglagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_C) #Nposlagtotal 
neglagaverage /= len(tracks_geo_region_C)

std_err_pos = np.std(all_ac_pos,axis=0,ddof=1)/np.sqrt(np.shape(all_ac_pos)[0])
std_err_neg = np.std(all_ac_neg,axis=0,ddof=1)/np.sqrt(np.shape(all_ac_neg)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
#plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err_pos,label='pos lag')
plt.errorbar(np.arange(0,min_track_length-4),neglagaverage[0:min_track_length-4],yerr=std_err_neg,label='neg lag')
plt.xlabel('lag (10 min)')
plt.title("Cross correlation velocity and polarity vectors {} for D/T {} - {}".format(region_name, np.round(bins[0],3), np.round(bins[1],3)))
plt.savefig('figures/binned_acf_figures/{}_vel_pol_crosscorr_avg_{}.png'.format(region, 'C'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_vel_pol_crosscorr_avg_{}.png \n'.format(region, 'C'))

#Autocorrelation speed
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region_A:
  track=df
  combined = pd.concat([track[['v']].reset_index(drop=True),track[['v']].reset_index(drop=True)], axis = 1 )
  combined = combined.dropna()
  #combined = combined.iloc[1::2]
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_A)# Nposlagtotal 

std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-1,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel('lag (5 min)')
plt.title("Autocorrelation speed {} for D/T {} - {}".format(region_name, np.round(bins[2],3), np.round(bins[3],3)))
plt.savefig('figures/binned_acf_figures/{}_speed_acf_avg_{}.png'.format(region, 'A'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_speed_acf_avg_{}.png \n'.format(region, 'A'))

#Autocorrelation speed
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region_B:
  track=df
  combined = pd.concat([track[['v']].reset_index(drop=True),track[['v']].reset_index(drop=True)], axis = 1 )
  combined = combined.dropna()
  #combined = combined.iloc[1::2]
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_B)# Nposlagtotal 

std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-1,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel('lag (5 min)')
plt.title("Autocorrelation speed {} for D/T {} - {}".format(region_name, np.round(bins[1],3), np.round(bins[2],3)))
plt.savefig('figures/binned_acf_figures/{}_speed_acf_avg_{}.png'.format(region, 'B'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_speed_acf_avg_{}.png \n'.format(region, 'B'))

#Autocorrelation speed
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region_C:
  track=df
  combined = pd.concat([track[['v']].reset_index(drop=True),track[['v']].reset_index(drop=True)], axis = 1 )
  combined = combined.dropna()
  #combined = combined.iloc[1::2]
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region_C)# Nposlagtotal 

std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-1,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel('lag (5 min)')
plt.title("Autocorrelation speed {} for D/T {} - {}".format(region_name, np.round(bins[0],3), np.round(bins[1],3)))
plt.savefig('figures/binned_acf_figures/{}_speed_acf_avg_{}.png'.format(region, 'C'))
plt.clf()
file_lines.append('figures/binned_acf_figures/{}_speed_acf_avg_{}.png \n'.format(region, 'C'))

#write lines to text file 
acf_fig_region.writelines(file_lines)
acf_fig_region.close() 