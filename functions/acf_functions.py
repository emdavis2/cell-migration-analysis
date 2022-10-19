import numpy as np

# for all acf functions, input is dfraw - a dataframe containing information for a single cell and the minimum length of the track specified in the analysis
# and returns acf for single cell up to min_track_length-4 (length is min_track_length-4 since calculating velocity involves smoothing window) normalized by first value of acf
# so all acf starts at 1 at lag 0


# xcorr_direction is used when calculating acf for velocity angle
def xcorr_direction(dfraw, min_track_length):
  df = dfraw.dropna()
  norms1 = (df.iloc[:,0]**2 + df.iloc[:,1]**2)**0.5 
  norms2 = (df.iloc[:,2]**2 + df.iloc[:,3]**2)**0.5
  v1x = np.asarray(df.iloc[:,0]/norms1)
  v1y = np.asarray(df.iloc[:,1]/norms1)
  v2x = np.asarray(df.iloc[:,2]/norms2)
  v2y = np.asarray(df.iloc[:,3]/norms2)

  length = len(df)
  poslagsmean=[]
  neglagsmean=[]
  Nposlags=[]
  Nneglags=[]
  for lag in range(length):
    poslags =  v2x[lag:length]*v1x[0:length-lag] + v2y[lag:length]*v1y[0:length-lag]
    neglags =  v2x[0:length-lag]*v1x[lag:length] + v2y[0:length-lag]*v1y[lag:length]
    poslagsmean.append(np.nanmean(poslags))
    neglagsmean.append(np.nanmean(neglags))
    Nposlags.append(sum(~np.isnan(poslags)))
    Nneglags.append(sum(~np.isnan(neglags)))

  return np.asarray(poslagsmean[0:min_track_length-4]), np.asarray(Nposlags[0:min_track_length-4]), np.asarray(neglagsmean[0:min_track_length-4]), np.asarray(Nneglags[0:min_track_length-4])


# xcorr_vector is used when calculating acf for polarity vector, polarity angle, and velocity
def xcorr_vector(dfraw, min_track_length):
  df = dfraw.dropna()
  v1x = np.asarray(df.iloc[:,0])
  v1y = np.asarray(df.iloc[:,1])
  v2x = np.asarray(df.iloc[:,2])
  v2y = np.asarray(df.iloc[:,3])

  length = len(df)
  poslagsmean=[]
  neglagsmean=[]
  Nposlags=[]
  Nneglags=[]
  for lag in range(length):
    poslags =  v2x[lag:length]*v1x[0:length-lag] + v2y[lag:length]*v1y[0:length-lag]
    neglags =  v2x[0:length-lag]*v1x[lag:length] + v2y[0:length-lag]*v1y[lag:length]
    poslagsmean.append(np.nanmean(poslags))
    neglagsmean.append(np.nanmean(neglags))
    Nposlags.append(sum(~np.isnan(poslags)))
    Nneglags.append(sum(~np.isnan(neglags)))

  return np.asarray(poslagsmean[0:min_track_length-4])/poslagsmean[0], np.asarray(Nposlags[0:min_track_length-4]), np.asarray(neglagsmean[0:min_track_length-4])/neglagsmean[0], np.asarray(Nneglags[0:min_track_length-4])


# xcorr is used when calculating abs-skew, speed, speed_x, and speed_y
def xcorr(dfraw, min_track_length):
  df = dfraw.dropna()
  v1 = np.asarray(df.iloc[:,0])  
  v2 = np.asarray(df.iloc[:,1])
  v1 = v1 - v1.mean()
  v2 = v2 - v2.mean()
  length = len(df)
  poslagsmean=[]
  neglagsmean=[]
  Nposlags=[]
  Nneglags=[]
  for lag in range(length):
    poslags =  v2[lag:]*v1[0:length-lag] 
    neglags =  v2[0:length-lag]*v1[lag:] 
    poslagsmean.append(np.nanmean(poslags))
    neglagsmean.append(np.nanmean(neglags))
    Nposlags.append(sum(~np.isnan(poslags)))
    Nneglags.append(sum(~np.isnan(neglags)))

  return np.asarray(poslagsmean[0:min_track_length-4])/poslagsmean[0], np.asarray(Nposlags[0:min_track_length-4]), np.asarray(neglagsmean[0:min_track_length-4])/neglagsmean[0], np.asarray(Nneglags[0:min_track_length-4])