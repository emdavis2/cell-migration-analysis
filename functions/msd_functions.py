import numpy as np

def calc_MSD(data,track_length):
  msd_allcells = []
  for cell in range(len(data)):
    length = len(data[cell])
    msd_onecell = []
    for lag in range(length):
      x = np.array(data[cell][:,0])
      y = np.array(data[cell][:,1])
      r = np.sqrt(x**2 + y**2)
      poslags =  (r[lag:] - r[0:length-lag])**2
      msd_onecell.append(np.average(poslags))
    msd_allcells.append(msd_onecell[0:track_length])
  return np.average(msd_allcells,axis=0)


def calc_MSD_sim(data,track_length):
  msd_allcells = []
  for cell in range(len(data)):
    length = len(data[cell])
    msd_onecell = []
    for lag in range(length):
      x = np.array(data[cell]['x'].to_list())
      y = np.array(data[cell]['y'].to_list())
      r = np.sqrt(x**2 + y**2)
      poslags =  (r[lag:] - r[0:length-lag])**2
      msd_onecell.append(np.average(poslags))
    msd_allcells.append(msd_onecell[0:track_length])
  return np.average(msd_allcells,axis=0)