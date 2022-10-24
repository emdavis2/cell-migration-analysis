from functions.acf_functions import *
from functions.compile_data_tracks_function import *

import ntpath
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

treatment = str(sys.argv[1])

min_track_length = int(sys.argv[2])

region = str(sys.argv[3])

tracks_region, tracks_geo_region, region_cells, region_endpointcells = compile_data_tracks(treatment, min_track_length, region)

#autocorrelation polarity vector
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  cospol = list( track[['abs-skew']].iloc[:,0] * np.cos(track[['polarity_angle']].iloc[:,0]) )
  sinpol = list( track[['abs-skew']].iloc[:,0] * np.sin(track[['polarity_angle']].iloc[:,0]) ) 
  normvect = pd.DataFrame( {'cospol': cospol , 'sinpol':sinpol })
  combined = pd.concat([normvect.reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean
poslagaverage /= len(tracks_geo_region) 

std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel("time lag")
plt.title("Autocorrelaton polarity_vector {}".format(region))
plt.savefig('figures/acf_figures/{}_polarity_vector_acf_avg.png'.format(region))
plt.clf()

#autocorrelation polarity angle
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  normvect = pd.DataFrame( {'cospolangle':list(np.cos(track[['polarity_angle']].iloc[:,0])) , 'sinpolangle':list(np.sin(track[['polarity_angle']].iloc[:,0])) })
  combined = pd.concat([normvect.reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean
poslagaverage /= len(tracks_geo_region) 

std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-0.6,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel("time lag")
plt.title("Autocorrelaton polarity_angle {}".format(region))
plt.savefig('figures/acf_figures/{}_polarity_anlge_acf_avg.png'.format(region))
plt.clf()

#Autocorrelation abs-skew
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  combined = pd.concat([track[['abs-skew']].reset_index(drop=True),track[['abs-skew']].reset_index(drop=True)], axis = 1 )
  combined = combined.dropna()
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)


  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region)  #Nposlagtotal 

std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-0.8,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel("time lag")
plt.title("Autocorrelation abs-skew {}".format(region))
plt.savefig('figures/acf_figures/{}_abs_skew_acf_avg.png'.format(region))
plt.clf()

#autocorrelation velocity angle
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_direction(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region) #Nposlagtotal 

std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-1,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel("time lag")
plt.title("Autocorrelaton velocity_angle {}".format(region))
plt.savefig('figures/acf_figures/{}_velocity_angle_acf_avg.png'.format(region))
plt.clf()

#autocorrelation velocity 
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region) #Nposlagtotal 

std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel("time lag")
plt.title("Autocorrelaton velocity {}".format(region))
plt.savefig('figures/acf_figures/{}_velocity_acf_avg.png'.format(region))
plt.clf()

#Autocorrelation speed
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  combined = pd.concat([track[['v']].reset_index(drop=True),track[['v']].reset_index(drop=True)], axis = 1 )
  combined = combined.dropna()
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region)# Nposlagtotal 

std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-1,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel("time lag")
plt.title("Autocorrelation speed {}".format(region))
plt.savefig('figures/acf_figures/{}_speed_acf_avg.png'.format(region))
plt.clf()


#Autocorrelation speed_x
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  combined = pd.concat([track[['vx']].reset_index(drop=True),track[['vx']].reset_index(drop=True)], axis = 1 )
  combined = combined.dropna()
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region)# Nposlagtotal 

std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-1,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel("time lag")
plt.title(" Autocorrelation speed_x {}".format(region))
plt.savefig('figures/acf_figures/{}_speed_x_acf_avg.png'.format(region))
plt.clf()

#Autocorrelation speed_y
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  combined = pd.concat([track[['vy']].reset_index(drop=True),track[['vy']].reset_index(drop=True)], axis = 1 )
  combined = combined.dropna()
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  all_ac.append(poslagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region)# Nposlagtotal 

std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-1,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel("time lag")
plt.title(" Autocorrelation speed_y {}".format(region))
plt.savefig('figures/acf_figures/{}_speed_y_acf_avg.png'.format(region))
plt.clf()