from functions.acf_functions import *
from functions.compile_data_tracks_function import *

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
for treat_name in tracks_geo_region:
  filtered_tracks_geo_region['firsthalf_{}'.format(treat_name)] = []
  filtered_tracks_geo_region['secondhalf_{}'.format(treat_name)] = []
  for df in tracks_geo_region[treat_name]:
    firsthalf_df = df.loc[df['frame']<108]
    secondhalf_df = df.loc[df['frame']>=108]
    if len(firsthalf_df) > 0:
      filtered_tracks_geo_region['firsthalf_{}'.format(treat_name)].append(firsthalf_df)
    if len(secondhalf_df) > 0:
      filtered_tracks_geo_region['secondhalf_{}'.format(treat_name)].append(secondhalf_df)


tracks_geo_region = filtered_tracks_geo_region


treatment_names = list(tracks_geo_region.keys())

#check to see if the path exists, if not make the directory
if not os.path.exists('sentinels'):
  os.mkdir('sentinels')

#clears out sentinel file if it exists
open('sentinels/ACF_figures.txt','w').close()
#create new sentinel file to write to
acf_fig_region = open('sentinels/ACF_figures.txt','w')
file_lines = []

#check to see if the path exists, if not make the directory
if not os.path.exists(save_path):
  os.mkdir(save_path)

#where figures are saved
figure_path = save_path + '/acf_figures'

#check to see if the path exists, if not make the directory
if not os.path.exists(figure_path):
  os.mkdir(figure_path)

for name in tracks_geo_region:

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
  '''
  #autocorrelation velocity angle
  poslagaverage = np.zeros(300)
  Nposlagtotal = np.zeros(300)
  all_ac = []
  for df in tracks_geo_region[name]:
    track=df
    combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
    poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_direction(combined, min_track_length)

    #remove nans here
    poslagsmean[np.isnan(poslagsmean)] = 0
    all_ac.append(poslagsmean)
    poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
    #Nposlagtotal[0:len(Nposlags)] += Nposlags
  poslagaverage /= len(tracks_geo_region[name]) #Nposlagtotal 

  std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

  #plt.plot(poslagaverage,label = "positive lag")
  plt.hlines(y=0,xmin=0,xmax=100,color='k')
  plt.xlim(0,min_track_length-4)
  plt.ylim(-1,1)
  plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
  plt.xlabel('lag (5 min)')
  plt.title("Autocorrelaton velocity_angle {}".format(name))
  plt.savefig('{}/{}_velocity_angle_acf_avg.png'.format(figure_path,name))
  plt.clf()
  file_lines.append('{}/{}_velocity_angle_acf_avg.png \n'.format(figure_path,name))
  '''
  #autocorrelation velocity 
  poslagaverage = np.zeros(300)
  Nposlagtotal = np.zeros(300)
  all_ac = []
  for df in tracks_geo_region[name]:
    track=df
    combined = pd.concat([track[['vx','vy']].reset_index(drop=True), track[['vx','vy']].reset_index(drop=True)], axis = 1 )
    #combined = combined.iloc[1::2]
    poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr_vector(combined, min_track_length)

    #remove nans here
    poslagsmean[np.isnan(poslagsmean)] = 0
    all_ac.append(poslagsmean)
    poslagaverage[0:len(poslagsmean)] += poslagsmean # Nposlags*poslagsmean
    #Nposlagtotal[0:len(Nposlags)] += Nposlags
  poslagaverage /= len(tracks_geo_region[name]) #Nposlagtotal 

  std_err = np.std(all_ac,axis=0,ddof=1)/np.sqrt(np.shape(all_ac)[0])

  plt.plot(poslagaverage,label = "positive lag")
  plt.hlines(y=0,xmin=0,xmax=100,color='k')
  plt.xlim(0,min_track_length-4)
  plt.ylim(-0.5,1)
  #plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
  plt.xlabel('lag (5 min)')
  plt.title("Autocorrelaton velocity {}".format(name))
  plt.savefig('{}/{}_velocity_acf_avg.png'.format(figure_path,name))
  plt.clf()
  file_lines.append('{}/{}_velocity_acf_avg.png \n'.format(figure_path,name))
  '''
  #Autocorrelation speed
  poslagaverage = np.zeros(300)
  Nposlagtotal = np.zeros(300)
  all_ac = []
  for df in tracks_geo_region[name]:
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
  poslagaverage /= len(tracks_geo_region[name])# Nposlagtotal 

  std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

  #plt.plot(poslagaverage,label = "positive lag")
  plt.hlines(y=0,xmin=0,xmax=100,color='k')
  plt.xlim(0,min_track_length-4)
  plt.ylim(-1,1)
  plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
  plt.xlabel('lag (5 min)')
  plt.title("Autocorrelation speed {}".format(name))
  plt.savefig('{}/{}_speed_acf_avg.png'.format(figure_path,name))
  plt.clf()
  file_lines.append('{}/{}_speed_acf_avg.png \n'.format(figure_path,name))


  #Autocorrelation speed_x
  poslagaverage = np.zeros(300)
  Nposlagtotal = np.zeros(300)
  all_ac = []
  for df in tracks_geo_region[name]:
    track=df
    combined = pd.concat([track[['vx']].reset_index(drop=True),track[['vx']].reset_index(drop=True)], axis = 1 )
    combined = combined.dropna()
    poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)

    #remove nans here
    poslagsmean[np.isnan(poslagsmean)] = 0
    all_ac.append(poslagsmean)
    poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
    #Nposlagtotal[0:len(Nposlags)] += Nposlags
  poslagaverage /= len(tracks_geo_region[name])# Nposlagtotal 

  std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

  #plt.plot(poslagaverage,label = "positive lag")
  plt.hlines(y=0,xmin=0,xmax=100,color='k')
  plt.xlim(0,min_track_length-4)
  plt.ylim(-1,1)
  plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
  plt.xlabel('lag (5 min)')
  plt.title(" Autocorrelation speed_x {}".format(name))
  plt.savefig('{}/{}_speed_x_acf_avg.png'.format(figure_path,name))
  plt.clf()
  file_lines.append('{}/{}_speed_x_acf_avg.png \n'.format(figure_path,name))

  #Autocorrelation speed_y
  poslagaverage = np.zeros(300)
  Nposlagtotal = np.zeros(300)
  all_ac = []
  for df in tracks_geo_region[name]:
    track=df
    combined = pd.concat([track[['vy']].reset_index(drop=True),track[['vy']].reset_index(drop=True)], axis = 1 )
    combined = combined.dropna()
    poslagsmean, Nposlags, neglagsmean, Nneglags = xcorr(combined, min_track_length)

    #remove nans here
    poslagsmean[np.isnan(poslagsmean)] = 0
    all_ac.append(poslagsmean)
    poslagaverage[0:len(poslagsmean)] += poslagsmean #Nposlags*poslagsmean
    #Nposlagtotal[0:len(Nposlags)] += Nposlags
  poslagaverage /= len(tracks_geo_region[name])# Nposlagtotal 

  std_err = np.std(all_ac,axis=0)/np.sqrt(np.shape(all_ac)[0])

  #plt.plot(poslagaverage,label = "positive lag")
  plt.hlines(y=0,xmin=0,xmax=100,color='k')
  plt.xlim(0,min_track_length-4)
  plt.ylim(-1,1)
  plt.errorbar(np.arange(0,min_track_length-4),poslagaverage[0:min_track_length-4],yerr=std_err)
  plt.xlabel('lag (5 min)')
  plt.title(" Autocorrelation speed_y {}".format(name))
  plt.savefig('{}/{}_speed_y_acf_avg.png'.format(figure_path,name))
  plt.clf()
  file_lines.append('{}/{}_speed_y_acf_avg.png \n'.format(figure_path,name))

  #cross correlation velocity vector and polarity vector
  poslagaverage = np.zeros(300)
  neglagaverage = np.zeros(300)
  #Nposlagtotal = np.zeros(300)
  all_ac_pos = []
  all_ac_neg = []
  for df in tracks_geo_region[name]:
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
  poslagaverage /= len(tracks_geo_region[name]) #Nposlagtotal 
  neglagaverage /= len(tracks_geo_region[name])

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
  plt.title("Cross correlation velocity and polarity vectors {}".format(name))
  plt.savefig('{}/{}_vel_pol_crosscorr_avg.png'.format(figure_path,name))
  plt.clf()
  file_lines.append('{}/{}_vel_pol_crosscorr_avg.png \n'.format(figure_path,name))

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

#write lines to text file 
acf_fig_region.writelines(file_lines)
acf_fig_region.close() 
'''