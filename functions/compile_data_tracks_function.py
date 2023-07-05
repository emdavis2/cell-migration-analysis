import ntpath
import os
import numpy as np
import pandas as pd

from functions.libraries.track_functions import *

# treatment path to folder containing pickle files of track data for region of interest (gel or glass)
# min_track_length is an integer of the minimum track length we are going to analyze
# region is a string of region of interest (stiff or glass)

#returns ind_track_regions (indicies of tracks of interest), tracks_region (list of x,y cell coordinates for each track), tracks_geo_region (list of dataframes for each track for region of interest with track parameters)
def compile_data_tracks(treatment, min_track_length, region, pixel_size):

    #Parameters to compute motion metrics
    center='approximate-medoid'
    #pixel_size=2*0.645 #um. factor of 2 because image was rescaled
    sampling_t=1 #sampling time (10min)

    #COMPILE TREATMENT TRACKS
    tracks =[]
    tracks_geo =[]

    for file in os.listdir(treatment):
        if file.endswith('.pkl'):
            tracksp,tracksgeop = read_tracks_aut([treatment+'/'+ file] ,pixel_size,center)
            if len(tracksp)==0: print("number of tracks read is 0")

            tracks += tracksp
            tracks_geo += tracksgeop

    lengths = [len(tracks_geo[i]) for i in range(len(tracks_geo)) ]
    print("number of tracks before filtering = ", len(tracks_geo))

    #remove tracks with less than 'mintracklength' points
    tracks = [track for track in tracks if len(track) >= min_track_length ]
    tracks_geo = [track for track in tracks_geo if len(track) >= min_track_length ]

    lengths = [len(tracks_geo[i]) for i in range(len(tracks_geo)) ]
    print("number of tracks after filtering = ", len(tracks_geo))

    #ADD FEATURES TO TRACKS

    #Add FMI, dFMI to tracks_geo
    for i in range(len(tracks_geo)):

        #get smooth tracks
        tracks_geo[i]['x_smooth3'] = pixel_size*tracks_geo[i][center+'x'].rolling(3,center=True).mean()
        #reverse sign of y due to inverted image axis
        tracks_geo[i]['y_smooth3'] = -pixel_size*tracks_geo[i][center+'y'].rolling(3,center=True).mean()
        tracks_geo[i]['dx_smooth3'] = tracks_geo[i]['x_smooth3'].diff()
        tracks_geo[i]['dy_smooth3'] = tracks_geo[i]['y_smooth3'].diff()
        #get velocity AT the timepoint
        tracks_geo[i]['vx'] = tracks_geo[i]['dx_smooth3'].rolling(2).mean().shift(-1)
        tracks_geo[i]['vy'] = tracks_geo[i]['dy_smooth3'].rolling(2).mean().shift(-1)
        tracks_geo[i]['v'] = ( tracks_geo[i]['vx']**2 + tracks_geo[i]['vy']**2 )**0.5

        tracks_geo[i]['velangle_smooth3']=np.arctan2(tracks_geo[i]['vy'], tracks_geo[i]['vx'] ) 
        
        tracks_geo[i]['polarity_turn'] =  tracks_geo[i]['polarity_angle'].diff()

        tracks_geo[i]['l_smooth3'] = ( tracks_geo[i]['dx_smooth3']**2 + tracks_geo[i]['dy_smooth3']**2 )**0.5

        #window quantities
        tracks_geo[i]['displ_win5'] = np.sqrt((tracks_geo[i]['dx_smooth3'].rolling(5,center=True).sum())**2 + (tracks_geo[i]['dy_smooth3'].rolling(5,center=True).sum())**2)
        tracks_geo[i]['pathlen_win5'] = tracks_geo[i]['l_smooth3'].rolling(5,center=True).sum()  
        tracks_geo[i]['DoTwin5'] = tracks_geo[i]['displ_win5']/tracks_geo[i]['pathlen_win5']
        tracks_geo[i]['abs-skew5'] = tracks_geo[i]['abs-skew'].rolling(5,center=True).mean()


        vec_step=tracks[i][1:]-tracks[i][0:-1]


        stepsize = np.sqrt((vec_step[:,0])**2 +(vec_step[:,1])**2)   

        #vector steps (velocity)
        #dvelocity
        dvel = vec_step[1:]-vec_step[0:-1]
        dvelmag = np.sqrt((dvel[:,0])**2 +(dvel[:,1])**2)    
        dvelang=np.arctan2(dvel[:,1],dvel[:,0])

        #Angle [-pi,pi]
        theta=np.arctan2(vec_step[:,1],vec_step[:,0])


        #cumulative path length, FMI, D/T (value from time zero to time j)   
        pathlength=np.zeros(len(tracks[i]) - 1)    
        fmi=np.zeros(len(tracks[i]) - 1)
        pmi=np.zeros(len(tracks[i]) - 1)
        DT=np.zeros(len(tracks[i]) - 1)        
        pathlength[0]=stepsize[0]
        #fmi[0] corresponds to the second timepoint
        fmi[0]=(tracks[i][1,0]- tracks[i][0,0])/pathlength[0]
        pmi[0]=(tracks[i][1,1]- tracks[i][0,1])/pathlength[0]
        DT[0]=( (tracks[i][1,0]-tracks[i][0,0])**2 + (tracks[i][1,1]-tracks[i][0,1])**2  )**0.5/pathlength[0]
        
        for j in range(1,len(tracks[i])-1): 
            pathlength[j] = pathlength[j-1] + stepsize[j]
            if pathlength[j] ==0:
                fmi[j]=0
                pmi[j]=0
                DT[j]=0
            else:
                #fmi[j]=tracks[i][j+1,0]
                DT[j]=( (tracks[i][j+1,0]-tracks[i][0,0])**2 + (tracks[i][j+1,1]-tracks[i][0,1])**2  )**0.5/pathlength[j]
                fmi[j]=(tracks[i][j+1,0] - tracks[i][0,0] )/pathlength[j]
                pmi[j]=(tracks[i][j+1,1] - tracks[i][0,1] )/pathlength[j]



        displmetrics = pd.DataFrame({'l':stepsize, 'theta':theta,'dvelmag':[np.nan]+list(dvelmag),
                                    'dvelang':[np.nan]+ list(dvelang), 'fmi':fmi, 'pmi':pmi, 'dovert':DT,
                                    'dx':vec_step[:,0], 'dy':vec_step[:,1]})   
        
        tracks_geo[i] = pd.concat([tracks_geo[i].reset_index(drop=True), displmetrics.reset_index(drop=True)], axis = 1 )

    #keep only tracks a particular region
    ind_tracks_region = [i for i in range(len(tracks_geo)) if tracks_geo[i]["region"][0] == region ]
    tracks_region = [tracks[i] for i in ind_tracks_region ]
    tracks_geo_region = [tracks_geo[i] for i in ind_tracks_region ]

    # Cell averaged metrics
    exclude_cols = pd.Series(['movie','frame','label','region'])
    colindices= ~tracks_geo[0].columns.isin(exclude_cols)

    #GET MEAN/STD GEOMETRY (for now pixel units)

    #make data frames with cell mean and std dev of not excluded shape metrics
    cellsshape = pd.DataFrame()
    cellsshape_std = pd.DataFrame()
    for i in range(len(tracks_geo)):
        #remove columns with duplicate names
        tracks_geo[i] = tracks_geo[i].loc[:,~tracks_geo[i].columns.duplicated()]
        #compute mean over cell
        colindices= ~tracks_geo[i].columns.isin(exclude_cols)
        means = tracks_geo[i].loc[:,colindices].mean().T
        #add gel-region
        means['region'] =tracks_geo[i]['region'].iloc[0]

        cellsshape = cellsshape.append(means, ignore_index=True)
        
        #get variation in shape features (STD)
        cellsshape_std = cellsshape_std.append(tracks_geo[i].loc[:,colindices].std().T, ignore_index=True)
    #set column names of std dev data frame
    cellsshape_std.columns = [feature+'_std' for feature in cellsshape_std.columns  ]

    cells = pd.concat([cellsshape.reset_index(drop=True), cellsshape_std.reset_index(drop=True)],1)

    #save as csv if wanted
    #cells.to_csv("cells.csv")

    #GET CELL-AVERAGED METRICS, get gel-region, experiment
    gelregion=[]
    experiment=[]
    trackid=[]
    absskew=[]
    mean_retr_norm_radii=[]

    for i in range(len(tracks_geo)):
        gelregion.append( tracks_geo[i]["region"].iloc[0])
        experiment.append(tracks_geo[i]["experiment"].iloc[0])
        trackid.append(tracks_geo[i]["track_id"].iloc[0])
        absskew.append(tracks_geo[i]["abs-skew"].mean())
        mean_retr_norm_radii.append(tracks_geo[i]["mean_retr_norm_radii"].mean())

    #GET TRACKS MOTION METRICS AND ADD GEL-REGION AND EXPERIMENT
    stepsizes,turns,meancoskturn,stderrcoskturn, tseries_stats, endpointcells = basic_stats(tracks,pixel_size,sampling_t)
    endpointcells['region']=gelregion
    endpointcells['experiment']=experiment
    endpointcells['track_id']=trackid
    endpointcells['abs-skew']=absskew
    endpointcells['mean_retr_norm_radii']=mean_retr_norm_radii

    #save as csv if wanted
    #endpointcells.to_csv("endpointtable.csv")

    return tracks_region, tracks_geo_region, cells, endpointcells