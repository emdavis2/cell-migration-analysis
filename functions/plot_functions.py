import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import f_oneway

#Get velocity for specified region using corresponding tracks_geo_region list of dataframes
#returns velocity (and x and y components) in terms of um/min as single list for entire region
def ExtractVelocity(tracks_geo_region, sampling_t):
    v_region = []
    vx_region = []
    vy_region = []
    for i in range(len(tracks_geo_region)):
        v_region.append(tracks_geo_region[i]['v'].dropna().tolist())
        vx_region.append(tracks_geo_region[i]['vx'].dropna().tolist())
        vy_region.append(tracks_geo_region[i]['vy'].dropna().tolist())

    v_region = (np.concatenate(v_region).ravel())/sampling_t
    vx_region = (np.concatenate(vx_region).ravel())/sampling_t
    vy_region = (np.concatenate(vy_region).ravel())/sampling_t

    return v_region, vx_region, vy_region

#Get abs-skew and solidity for specified region using corresponding tracks_geo_region list of dataframes
#returns absskew and solidity as single lists for entire region
def ExtractAbsskewSolidity(tracks_geo_region):
    absskew_region = []
    solidity_region = []
    for i in range(len(tracks_geo_region)):
        absskew_region.append(tracks_geo_region[i]['abs-skew'].dropna().tolist())
        solidity_region.append(tracks_geo_region[i]['solidity'].dropna().tolist())

    absskew_region = np.concatenate(absskew_region).ravel()
    solidity_region = np.concatenate(solidity_region).ravel()

    return absskew_region, solidity_region


#Plot a histogram of data of interest for one region
def HistogramPlot(data_to_plot, region, metric, bin_num, save_path, sentinel_name):
    plt.hist(data_to_plot,bins=bin_num)
    plt.title('{} for {}'.format(metric, region))
    plt.xlabel('{}'.format(metric))
    plt.ylabel('counts')
    plt.savefig('{}/{}_hist_{}.png'.format(save_path, metric, region))
    plt.clf()
    sentinel_name.append('{}/{}_hist_{}.png \n'.format(save_path, metric, region))


#Plot boxplots comparing all regions for specified metric
def BoxplotPlot(region_list, metric_list, metric, save_path, sentinel_path, pval_file_name):
    pval_file_name.append('\n {} \n'.format(metric))
    data_bp = {}
    for index,region in enumerate(region_list):
        data_bp[region] = metric_list[region]
        tstat, pval =  f_oneway(metric_list[region_list[index]],metric_list[region_list[index-1]])
        pval_file_name.append('{} and {}: {} \n'.format(region_list[index], region_list[index-1], pval))
    data_boxplot = pd.DataFrame({ key:pd.Series(value) for key, value in data_bp.items() })
    sns.boxplot(data=data_boxplot)
    plt.xlabel("Source")
    plt.ylabel("{}".format(metric))
    plt.savefig('{}/{}_boxplot.png'.format(save_path, metric))
    plt.clf()
    sentinel_path.append('{}/{}_boxplot.png \n'.format(save_path, metric))

def PlotTrackLengthHist(tracks_geo_region, save_path, sentinel_name):
    for region in tracks_geo_region:
        lengths = [len(tracks_geo_region[region][i]) for i in range(len(tracks_geo_region[region])) ]
        num_tracks = len(tracks_geo_region[region])
        plt.hist(lengths)
        plt.xlabel('track length (number of frames)')
        plt.ylabel('counts')
        plt.title('Track lengths for {} data ({} tracks)'.format(region, num_tracks))
        plt.savefig('{}/track_lengths_data_{}.png'.format(save_path, region))
        plt.clf()
        sentinel_name.append('{}/track_lengths_data_{}.png'.format(save_path, region))


def PlotTrajectories(tracks_region, save_path, sentinel_name):
    for region in tracks_region:
        for df in tracks_region[region]:
            plt.plot(df[:,0],df[:,1])
        plt.xlabel(('x position ($\mu$m)'))
        plt.ylabel(('y position ($\mu$m)'))
        plt.title('Trajectories for {} Data'.format(region))
        plt.savefig('{}/trajectories_data_{}.png'.format(save_path, region))
        plt.clf()
        sentinel_name.append('{}/trajectories_data_{}.png \n'.format(save_path, region))
