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

tracks_region, tracks_geo_region, region_cells, region_endpointcells = compile_data_tracks(treatment, min_track_length, region, pixel_size)

if region == 'stiff':
  region_name = 'gel'
else:
  region_name = region

#clears out sentinel file if it exists
open('sentinels/ACF_figures_{}.txt'.format(region),'w').close()
#create new sentinel file to write to
acf_fig_region = open('sentinels/ACF_figures_{}.txt'.format(region),'w')
file_lines = []

# #autocorrelation polarity vector
# poslagaverage = np.zeros(300)
# Nposlagtotal = np.zeros(300)
# all_ac = []
# for df in tracks_geo_region:
#   track=df
#   cospol = list( track[['abs-skew']].iloc[:,0] * np.cos(track[['polarity_angle']].iloc[:,0]) )
#   sinpol = list( track[['abs-skew']].iloc[:,0] * np.sin(track[['polarity_angle']].iloc[:,0]) ) 
#   normvect = pd.DataFrame( {'cospol': cospol , 'sinpol':sinpol })
#   combined = pd.concat([normvect.reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
#   poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

#   #remove nans here
#   poslagsmean[np.isnan(poslagsmean)] = 0
#   all_ac.append(poslagsmean)
#   poslagaverage[0:len(poslagsmean)] += poslagsmean
# poslagaverage /= len(tracks_geo_region) 

# std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

# #plt.plot(poslagaverage,label = "positive lag")
# plt.hlines(y=0,xmin=0,xmax=100,color='k')
# plt.xlim(0,min_track_length-4)
# plt.ylim(-0.5,1)
# plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
# plt.xlabel('lag (10 min)')
# plt.title("Autocorrelaton polarity_vector {}".format(region_name))
# plt.savefig('figures/acf_figures/{}_polarity_vector_acf_avg.png'.format(region))
# plt.clf()
# file_lines.append('figures/acf_figures/{}_polarity_vector_acf_avg.png \n'.format(region))

# #autocorrelation polarity angle
# poslagaverage = np.zeros(300)
# Nposlagtotal = np.zeros(300)
# all_ac = []
# for df in tracks_geo_region:
#   track=df
#   normvect = pd.DataFrame( {'cospolangle':list(np.cos(track[['polarity_angle']].iloc[:,0])) , 'sinpolangle':list(np.sin(track[['polarity_angle']].iloc[:,0])) })
#   combined = pd.concat([normvect.reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
#   poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

#   #remove nans here
#   poslagsmean[np.isnan(poslagsmean)] = 0
#   all_ac.append(poslagsmean)
#   poslagaverage[0:len(poslagsmean)] += poslagsmean
# poslagaverage /= len(tracks_geo_region) 

# std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

# #plt.plot(poslagaverage,label = "positive lag")
# plt.hlines(y=0,xmin=0,xmax=100,color='k')
# plt.xlim(0,min_track_length-4)
# plt.ylim(-0.6,1)
# plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
# plt.xlabel('lag (5 min)')
# plt.title("Autocorrelaton polarity_angle {}".format(region_name))
# plt.savefig('figures/acf_figures/{}_polarity_angle_acf_avg.png'.format(region))
# plt.clf()
# file_lines.append('figures/acf_figures/{}_polarity_angle_acf_avg.png \n'.format(region))

# #Autocorrelation abs-skew
# poslagaverage = np.zeros(300)
# Nposlagtotal = np.zeros(300)
# all_ac = []
# for df in tracks_geo_region:
#   track=df
#   combined = pd.concat([track[['abs-skew']].reset_index(drop=True),track[['abs-skew']].reset_index(drop=True)], axis = 1 )
#   combined = combined.dropna()
#   poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)


#   #remove nans here
#   poslagsmean[np.isnan(poslagsmean)] = 0
#   all_ac.append(poslagsmean)
#   poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
#   #Nposlagtotal[0:len(Nposlags)] += Nposlags
# poslagaverage /= len(tracks_geo_region)  #Nposlagtotal 

# std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

# #plt.plot(poslagaverage,label = "positive lag")
# plt.hlines(y=0,xmin=0,xmax=100,color='k')
# plt.xlim(0,min_track_length-4)
# plt.ylim(-0.8,1)
# plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
# plt.xlabel('lag (5 min)')
# plt.title("Autocorrelation abs-skew {}".format(region_name))
# plt.savefig('figures/acf_figures/{}_abs_skew_acf_avg.png'.format(region))
# plt.clf()
# file_lines.append('figures/acf_figures/{}_abs_skew_acf_avg.png \n'.format(region))

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
plt.xlabel('lag (5 min)')
plt.title("Autocorrelaton velocity_angle {}".format(region_name))
plt.savefig('figures/acf_figures/{}_velocity_angle_acf_avg.png'.format(region))
plt.clf()
file_lines.append('figures/acf_figures/{}_velocity_angle_acf_avg.png \n'.format(region))

#autocorrelation velocity 
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
  track=df
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
  #combined = combined.iloc[1::2]
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
plt.xlabel('lag (5 min)')
plt.title("Autocorrelaton velocity {}".format(region_name))
plt.savefig('figures/acf_figures/{}_velocity_acf_avg.png'.format(region))
plt.clf()
file_lines.append('figures/acf_figures/{}_velocity_acf_avg.png \n'.format(region))

#Autocorrelation speed
poslagaverage = np.zeros(300)
Nposlagtotal = np.zeros(300)
all_ac = []
for df in tracks_geo_region:
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
poslagaverage /= len(tracks_geo_region)# Nposlagtotal 

std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
plt.ylim(-1,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
plt.xlabel('lag (5 min)')
plt.title("Autocorrelation speed {}".format(region_name))
plt.savefig('figures/acf_figures/{}_speed_acf_avg.png'.format(region))
plt.clf()
file_lines.append('figures/acf_figures/{}_speed_acf_avg.png \n'.format(region))


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
plt.xlabel('lag (5 min)')
plt.title(" Autocorrelation speed_x {}".format(region_name))
plt.savefig('figures/acf_figures/{}_speed_x_acf_avg.png'.format(region))
plt.clf()
file_lines.append('figures/acf_figures/{}_speed_x_acf_avg.png \n'.format(region))

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
plt.xlabel('lag (5 min)')
plt.title(" Autocorrelation speed_y {}".format(region_name))
plt.savefig('figures/acf_figures/{}_speed_y_acf_avg.png'.format(region))
plt.clf()
file_lines.append('figures/acf_figures/{}_speed_y_acf_avg.png \n'.format(region))

#write lines to text file 
acf_fig_region.writelines(file_lines)
acf_fig_region.close() 

#cross correlation velocity vector and polarity vector
poslagaverage = np.zeros(300)
neglagaverage = np.zeros(300)
#Nposlagtotal = np.zeros(300)
all_ac_pos = []
all_ac_neg = []
for df in tracks_geo_region:
  track=df
  cospol = list( np.cos(track[['polarity_angle']].iloc[:,0]) )
  sinpol = list( np.sin(track[['polarity_angle']].iloc[:,0]) ) 
  normvect = pd.DataFrame( {'cospol': cospol , 'sinpol':sinpol })
  combined = pd.concat([track[['vx','vy']].reset_index(drop=True), normvect.reset_index(drop=True)], axis = 1 )
  poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_direction(combined, min_track_length)

  #remove nans here
  poslagsmean[np.isnan(poslagsmean)] = 0
  neglagsmean[np.isnan(neglagsmean)] = 0
  all_ac_pos.append(poslagsmean)
  all_ac_neg.append(neglagsmean)
  poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
  neglagaverage[0:len(neglagsmean)] += neglagsmean
  #Nposlagtotal[0:len(Nposlags)] += Nposlags
poslagaverage /= len(tracks_geo_region) #Nposlagtotal 
neglagaverage /= len(tracks_geo_region)

std_err_pos = np.std(all_ac_pos,axis=0,ddof=1)/np.sqrt(np.shape(all_ac_pos)[0])
std_err_neg = np.std(all_ac_neg,axis=0,ddof=1)/np.sqrt(np.shape(all_ac_neg)[0])

#plt.plot(poslagaverage,label = "positive lag")
plt.hlines(y=0,xmin=0,xmax=100,color='k')
plt.xlim(0,min_track_length-4)
#plt.ylim(-0.5,1)
plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err_pos,label='pos lag')
plt.errorbar(np.arange(0,min_track_length-4),neglagaverage[0:min_track_length-4],yerr=std_err_neg,label='neg lag')
plt.xlabel('lag (10 min)')
plt.legend()
plt.title("Cross correlation velocity and polarity vectors {}".format(region_name))
plt.savefig('figures/acf_figures/{}_vel_pol_crosscorr_avg.png'.format(region))
plt.clf()
file_lines.append('figures/acf_figures/{}_vel_pol_crosscorr_avg.png \n'.format(region))

# #cross correlation polarity vector and velocity vector
# poslagaverage = np.zeros(300)
# Nposlagtotal = np.zeros(300)
# all_ac = []
# for df in tracks_geo_region:
#   track=df
#   cospol = list( track[['abs-skew']].iloc[:,0] * np.cos(track[['polarity_angle']].iloc[:,0]) )
#   sinpol = list( track[['abs-skew']].iloc[:,0] * np.sin(track[['polarity_angle']].iloc[:,0]) ) 
#   normvect = pd.DataFrame( {'cospol': cospol , 'sinpol':sinpol })
#   combined = pd.concat([normvect.reset_index(drop=True),track[['vx','vy']].reset_index(drop=True)], axis = 1 )
#   poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

#   #remove nans here
#   poslagsmean[np.isnan(poslagsmean)] = 0
#   all_ac.append(poslagsmean)
#   poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
#   #Nposlagtotal[0:len(Nposlags)] += Nposlags
# poslagaverage /= len(tracks_geo_region) #Nposlagtotal 

# std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

# #plt.plot(poslagaverage,label = "positive lag")
# plt.hlines(y=0,xmin=0,xmax=100,color='k')
# plt.xlim(0,min_track_length-4)
# #plt.ylim(-0.5,1)
# plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
# plt.xlabel('lag (10 min)')
# plt.title("Cross correlation polarity and velocity vectors {}".format(region_name))
# plt.savefig('figures/acf_figures/{}_pol_vel_crosscorr_avg.png'.format(region))
# plt.clf()
# file_lines.append('figures/acf_figures/{}_pol_vel_crosscorr_avg.png \n'.format(region))