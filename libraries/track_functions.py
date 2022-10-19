import numpy as np
import random
import scipy.io
import os.path
from os import listdir
from os.path import isfile, join
import re
import csv
from copy import deepcopy
from scipy.stats import ks_2samp
import pandas as pd
from numpy import linalg as LA
import pickle

def read_tracks(path,files,pixel_size):
    tracks=[]
    for filename in files:
        csvdata=csv.reader(open(path+filename))
        data=["\t".join(i)for i in csvdata]
        raw_tracks=np.genfromtxt(data, delimiter="\t",skip_header=1)
        current_track=[]
        j=raw_tracks[0,1] #first track number
        for i in range(len(raw_tracks)):
            if raw_tracks[i,1] != j:
                #reset tracks index
                j=raw_tracks[i,1]
                #append current track to list of tracks
                tracks.append(current_track)
                #start a new track
                current_track=[pixel_size*raw_tracks[i,3:5]]
            else:
                #add coordinates to current track
                current_track.append(pixel_size*raw_tracks[i,3:5])
        #add last track to tracks
        tracks.append(current_track)

    #set initial point at the origin
    nptracks=[]
    #print (str(len(tracks))+' tracks')
    for i in range(len(tracks)):
        nptracks.append(np.asarray(tracks[i]) - tracks[i][0]   )

    tracks=nptracks

    #hapto, duro cells move upwards in movies,
    #but because y axis increases downwards, measured angle is as moving downwards (-pi/2)
    #therofore rotate data 90 degrees counterclock wise so that gradient is towards x axis (0 degrees)
    for i in range(len(tracks)):
        tempx = deepcopy(tracks[i][:,0])
        tracks[i][:,0]= -tracks[i][:,1]
        tracks[i][:,1]= tempx

    return tracks

def read_tracks_aut(files,pixel_size,center):
    tracks=[]
    tracksgeo=[]
    for filename in files:

        with open(filename, 'rb') as handle:
            tracksfile = pickle.load(handle, encoding='latin1')
        #print('number of tracks',len(tracksfile))
        for i in range(len(tracksfile)):
            #print('track length',len(tracksfile[i]))
            #print(tracksfile[i])
            coord=np.vstack((tracksfile[i][center+'x'],tracksfile[i][center+'y']))
            #tracksgeo.append(np.asarray(tracksfile[i]['area']))
            tracksgeo.append(tracksfile[i])
            #coor shape 2,ntimes
            tracks.append(np.transpose(coord)*pixel_size)
    #print (str(len(tracks))+' tracks')

    #set initial point at the origin and convert to np.array
    nptracks=[]
    #print (str(len(tracks))+' tracks')
    for i in range(len(tracks)):
        nptracks.append(np.asarray(tracks[i]) - tracks[i][0]   )

    tracks=nptracks
    #hapto, duro cells move upwards in movies,
    #but because y axis increases downwards I record them as moving downwards (-pi/2)
    #therofore rotate 90 degrees counterclock wise so that gradient is toward x axis
    for i in range(len(tracks)):
        tempx = deepcopy(tracks[i][:,0])
        tracks[i][:,0]= -tracks[i][:,1]
        tracks[i][:,1]= tempx

    return tracks, tracksgeo

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def smooth_tracks(tracks,n):
    newtracks=np.copy(tracks)
    for i in range(len(tracks)):
        if len(newtracks[i]) > n:
            newtracks[i][:,0]=smooth(newtracks[i][:,0],n)
            newtracks[i][:,1]=smooth(newtracks[i][:,1],n)
            newtracks[i]=newtracks[i][0:-2]
    return newtracks

def basic_stats(tracks,pixel_size,sampling_t):
    #input tracks must be centered at origin and with distance units of um
    stepsizes=np.empty(0)
    steps0=[]
    stepsizesturn=np.empty(0)
    turns=np.empty(0)
    thetas=np.empty(0)
    maxkturn=40
    kturns=[]
    regstats=[]
    for k in range(maxkturn):
        kturns.append(np.empty(0))

    siz=200
    sum_pathlengths=np.zeros(siz)
    sum_squared_pathlengths=np.zeros(siz)
    sum_fmis=np.zeros(siz)
    sum_squared_fmis=np.zeros(siz)
    sum_pmis=np.zeros(siz)
    sum_squared_pmis=np.zeros(siz)
    sum_DTs=np.zeros(siz)
    sum_squared_DTs=np.zeros(siz)
    nlengths=np.zeros(siz)

    for i in range(len(tracks)):
    #    print len(tracks[i])
        stepsize = np.sqrt((tracks[i][0:-1,0]-tracks[i][1:,0])**2 +(tracks[i][0:-1,1]-tracks[i][1:,1])**2)
        stepsizes=np.concatenate((stepsizes,stepsize))
        stepsizesturn=np.concatenate((stepsizesturn,stepsize[1:]))

        steps0.append(np.sum(stepsize==0))

        #vector steps
        vec_step=tracks[i][1:]-tracks[i][0:-1]

        #Angle [-pi,pi]
        theta=np.arctan2(vec_step[:,1],vec_step[:,0])
        #set nan when the step size is too small
        theta[(stepsize <= np.sqrt(2)*pixel_size )] = np.nan
        #theta[(stepsize <= pixel_size )] = np.nan
        thetas=np.concatenate((thetas,theta))

        #turning angle
        turnp=theta[1:]-theta[0:-1]
        #turning angle [-pi,pi]
        turn=np.arctan2(np.sin(turnp),np.cos(turnp))
        turns=np.concatenate((turns,turn))

        #turning angle after k steps
        for k in range(1,maxkturn+1):
            turnkp=theta[k:]-theta[0:-k]
            turnk=np.arctan2(np.sin(turnkp),np.cos(turnkp))
            kturns[k-1]=np.concatenate((kturns[k-1],turnk))

        #path length, FMI, D/T
        pathlength=np.zeros(len(tracks[i]))
        fmi=np.zeros(len(tracks[i]))
        pmi=np.zeros(len(tracks[i]))
        DT=np.zeros(len(tracks[i]))
        #pathlength[0]=stepsize[0]
        pathlength[0]=0
        fmi[0]=0
        pmi[0]=0
        DT[0]=0
        for j in range(1,len(tracks[i])):
            pathlength[j] = pathlength[j-1] + stepsize[j-1]
            if pathlength[j] ==0:
                fmi[j]=0
                pmi[j]=0
                DT[j]=0
            else:
                #fmi[j]=tracks[i][j+1,0]
                DT[j]=((tracks[i][j,0] - tracks[i][0,0])**2 + (tracks[i][j,1]-tracks[i][0,1])**2  )**0.5/pathlength[j]
                fmi[j]=(tracks[i][j,0] - tracks[i][0,0])/pathlength[j]
                pmi[j]=(tracks[i][j,1] - tracks[i][0,1])/pathlength[j]

         #regular statistics with last element:
        T = pathlength[-1]
        D = LA.norm(tracks[i][-1,0:2] - tracks[i][0,0:2] )
        yfinal = tracks[i][-1,1] - tracks[i][0,1]
        xfinal = tracks[i][-1,0] - tracks[i][0,0]
        angle = np.arctan2(yfinal,xfinal)
        #regstats: len, T, speed(um/hr), D, D/T, FMI, PMI, angle
        regstats.append([len(tracks[i]), T, T/sampling_t/(len(tracks[i])-1),
                         D, D/T, xfinal/T, yfinal/T , angle])

        nlengths[0:len(pathlength)] += np.ones(len(pathlength))
        sum_pathlengths[0:len(pathlength)] += pathlength
        sum_squared_pathlengths[0:len(pathlength)] += pathlength**2
        sum_fmis[0:len(pathlength)] += fmi
        sum_squared_fmis[0:len(pathlength)] += fmi**2
        sum_pmis[0:len(pathlength)] += pmi
        sum_squared_pmis[0:len(pathlength)] += pmi**2
        sum_DTs[0:len(pathlength)] += DT
        sum_squared_DTs[0:len(pathlength)] += DT**2

    endpointstatspercell = pd.DataFrame(data=np.asarray(regstats),
                                        index=range(len(regstats)),
                                        columns=['length','T','speed','D',
                                                 'DoverT','FMI','PMI','angle'])

    meanpathlength=sum_pathlengths/nlengths
    stddevpathlength = np.sqrt(sum_squared_pathlengths/nlengths - meanpathlength**2 )
    stderrpathlength = stddevpathlength/nlengths**0.5
    meanfmi=sum_fmis/nlengths
    stddevfmi = np.sqrt(sum_squared_fmis/nlengths - meanfmi**2 )
    stderrfmi = stddevfmi/nlengths**0.5
    meanpmi=sum_pmis/nlengths
    stddevpmi = np.sqrt(sum_squared_pmis/nlengths - meanpmi**2 )
    stderrpmi = stddevpmi/nlengths**0.5
    meanDT=sum_DTs/nlengths
    stddevDT = np.sqrt(sum_squared_DTs/nlengths - meanDT**2 )
    stderrDT = stddevDT/nlengths**0.5

    meancoskturn=[1]
    stderrcoskturn=[0]
    for i in range(len(kturns)):
        meancoskturn.append(np.nanmean(np.cos(kturns[i])))
        stderrcoskturn.append(np.nanstd(np.cos(kturns[i]))/np.sum(~np.isnan(kturns[i]))**0.5 )
        #print len(kturns[i])
    stderrcoskturn=np.asarray(stderrcoskturn)
    meancoskturn = np.asarray(meancoskturn)

    tseries_stats=pd.DataFrame({'Tvst':meanpathlength,'stderrTvst':stderrpathlength,
                                 'FMIvst':meanfmi, 'stderrFMIvst':stderrfmi,
                                 'PMIvst':meanpmi, 'stderrPMIvst': stderrpmi,
                                 'DoverTvst':meanDT, 'stderrDoverTvst':stderrDT })

    return stepsizes,turns,meancoskturn,stderrcoskturn, tseries_stats, endpointstatspercell

def basic_stats_window(tracks,pixel_size,sampling_t):
    #input tracks must be centered at origin and with distance units of um
    stepsizes=np.empty(0)
    steps0=[]
    stepsizesturn=np.empty(0)
    turns=np.empty(0)
    thetas=np.empty(0)
    maxkturn=40
    kturns=[]
    regstats=[]
    motion_metrics={}
    for k in range(maxkturn):
        kturns.append(np.empty(0))

    siz=200
    sum_pathlengths=np.zeros(siz)
    sum_squared_pathlengths=np.zeros(siz)
    sum_fmis=np.zeros(siz)
    sum_squared_fmis=np.zeros(siz)
    sum_pmis=np.zeros(siz)
    sum_squared_pmis=np.zeros(siz)
    sum_DTs=np.zeros(siz)
    sum_squared_DTs=np.zeros(siz)
    nlengths=np.zeros(siz)

    Pathlength=[]
    FMI=[]
    PMI=[]
    DDT=[]
    Speed=[]

    for i in range(len(tracks)):
    #    print len(tracks[i])
        stepsize = np.sqrt((tracks[i][0:-1,0]-tracks[i][1:,0])**2 +(tracks[i][0:-1,1]-tracks[i][1:,1])**2)
        stepsizes=np.concatenate((stepsizes,stepsize))
        stepsizesturn=np.concatenate((stepsizesturn,stepsize[1:]))

        steps0.append(np.sum(stepsize==0))

        #vector steps
        vec_step=tracks[i][1:]-tracks[i][0:-1]

        #Angle [-pi,pi]
        theta=np.arctan2(vec_step[:,1],vec_step[:,0])
        #set nan when the step size is too small
        theta[(stepsize <= np.sqrt(2)*pixel_size )] = np.nan
        #theta[(stepsize <= pixel_size )] = np.nan
        thetas=np.concatenate((thetas,theta))

        #turning angle
        turnp=theta[1:]-theta[0:-1]
        #turning angle [-pi,pi]
        turn=np.arctan2(np.sin(turnp),np.cos(turnp))
        turns=np.concatenate((turns,turn))

        #turning angle after k steps
        for k in range(1,maxkturn+1):
            turnkp=theta[k:]-theta[0:-k]
            turnk=np.arctan2(np.sin(turnkp),np.cos(turnkp))
            kturns[k-1]=np.concatenate((kturns[k-1],turnk))

        #path length, FMI, D/T
        pathlength=np.zeros(len(tracks[i]))
        fmi=np.zeros(len(tracks[i]))
        pmi=np.zeros(len(tracks[i]))
        DT=np.zeros(len(tracks[i]))
        speed=np.zeros(len(tracks[i]))
        #pathlength[0]=stepsize[0]
        pathlength[0]=0
        fmi[0]=0
        pmi[0]=0
        DT[0]=0
        speed[0]=0
        for j in range(1,len(tracks[i])):
            pathlength[j] = pathlength[j-1] + stepsize[j-1]
            if pathlength[j] ==0:
                fmi[j]=0
                pmi[j]=0
                DT[j]=0
                speed[j]=0
            else:
                #fmi[j]=tracks[i][j+1,0]
                DT[j]=((tracks[i][j,0] - tracks[i][0,0])**2 + (tracks[i][j,1]-tracks[i][0,1])**2  )**0.5/pathlength[j]
                fmi[j]=(tracks[i][j,0] - tracks[i][0,0])/pathlength[j]
                pmi[j]=(tracks[i][j,1] - tracks[i][0,1])/pathlength[j]
                speed[j]=(pathlength[j]/sampling_t/j)

        DDT.append(DT)
        FMI.append(fmi)
        PMI.append(pmi)
        Speed.append(speed)

         #regular statistics with last element:
        T = pathlength[-1]
        D = LA.norm(tracks[i][-1,0:2] - tracks[i][0,0:2] )
        yfinal = tracks[i][-1,1] - tracks[i][0,1]
        xfinal = tracks[i][-1,0] - tracks[i][0,0]
        angle = np.arctan2(yfinal,xfinal)
        #regstats: len, T, speed(um/hr), D, D/T, FMI, PMI, angle
        regstats.append([len(tracks[i]), T, T/sampling_t/(len(tracks[i])-1),
                         D, D/T, xfinal/T, yfinal/T , angle])

        #motion_metrics
        motion_metrics['DT']=DDT
        motion_metrics['FMI']=FMI
        motion_metrics['PMI']=PMI
        motion_metrics['speed']=Speed

        nlengths[0:len(pathlength)] += np.ones(len(pathlength))
        sum_pathlengths[0:len(pathlength)] += pathlength
        sum_squared_pathlengths[0:len(pathlength)] += pathlength**2
        sum_fmis[0:len(pathlength)] += fmi
        sum_squared_fmis[0:len(pathlength)] += fmi**2
        sum_pmis[0:len(pathlength)] += pmi
        sum_squared_pmis[0:len(pathlength)] += pmi**2
        sum_DTs[0:len(pathlength)] += DT
        sum_squared_DTs[0:len(pathlength)] += DT**2

    endpointstatspercell = pd.DataFrame(data=np.asarray(regstats),
                                        index=range(len(regstats)),
                                        columns=['length','T','speed','D',
                                                 'DoverT','FMI','PMI','angle'])

    meanpathlength=sum_pathlengths/nlengths
    stddevpathlength = np.sqrt(sum_squared_pathlengths/nlengths - meanpathlength**2 )
    stderrpathlength = stddevpathlength/nlengths**0.5
    meanfmi=sum_fmis/nlengths
    stddevfmi = np.sqrt(sum_squared_fmis/nlengths - meanfmi**2 )
    stderrfmi = stddevfmi/nlengths**0.5
    meanpmi=sum_pmis/nlengths
    stddevpmi = np.sqrt(sum_squared_pmis/nlengths - meanpmi**2 )
    stderrpmi = stddevpmi/nlengths**0.5
    meanDT=sum_DTs/nlengths
    stddevDT = np.sqrt(sum_squared_DTs/nlengths - meanDT**2 )
    stderrDT = stddevDT/nlengths**0.5

    meancoskturn=[1]
    stderrcoskturn=[0]
    for i in range(len(kturns)):
        meancoskturn.append(np.nanmean(np.cos(kturns[i])))
        stderrcoskturn.append(np.nanstd(np.cos(kturns[i]))/np.sum(~np.isnan(kturns[i]))**0.5 )
        #print len(kturns[i])
    stderrcoskturn=np.asarray(stderrcoskturn)
    meancoskturn = np.asarray(meancoskturn)

    tseries_stats=pd.DataFrame({'Tvst':meanpathlength,'stderrTvst':stderrpathlength,
                                 'FMIvst':meanfmi, 'stderrFMIvst':stderrfmi,
                                 'PMIvst':meanpmi, 'stderrPMIvst': stderrpmi,
                                 'DoverTvst':meanDT, 'stderrDoverTvst':stderrDT })

    return stepsizes,turns,meancoskturn,stderrcoskturn, tseries_stats, endpointstatspercell, motion_metrics

def bprw_2angles_sim(p,stepsizesdata):
    ntracks=1000
    nsteps=6*3 + 41 #accounting for 3 hours previous to recording (sims sampled every 10 mins)
    nturns=nsteps*ntracks
    kappamove,kappapol,w1,tx = p
    mumove=0
    mupol=0
    #grad=np.pi*3/2 #angle of orientation of the gradient
    grad=0
    #random.seed(1)
    #GENERATE ANGLES FROM PROBABILITY DENSITY FUNCTION TO USE IN SIMULATIONS
    noisepol=np.random.vonmises(mupol, kappapol,nturns)
    noisemove=np.random.vonmises(mumove, kappamove, nturns)
    #PATH SIMULATIONS
    tracks=[]
    counter=0
    for j in range(ntracks):
        #initial step
        thetapol=[np.random.uniform(-np.pi,np.pi)]
        thetasim=[thetapol[0]]
        stepsizesim=[random.choice(stepsizesdata)];
        current_track=[[0,0]]
        for i in range(1,nsteps+1):
            sum_grad_prevpol=np.arctan2((1-tx)*np.sin(thetapol[i-1]) + tx*np.sin(grad) , (1-tx)*np.cos(thetapol[i-1]) + tx*np.cos(grad) )
            #get theta pol i
            thetapol.append(sum_grad_prevpol + noisepol[counter])
            #short time scale process
            sum_prevth_thpoli=np.arctan2((1-w1)*np.sin(thetasim[i-1]) + w1*np.sin(thetapol[i-1]) , (1-w1)*np.cos(thetasim[i-1]) + w1*np.cos(thetapol[i-1]))
            #get theta i
            thetasim.append(sum_prevth_thpoli + noisemove[counter])
            #get step size i
            stepsizesim.append(random.choice(stepsizesdata))
            current_track.append([current_track[i-1][0] + stepsizesim[i]*np.cos(thetasim[i]) , current_track[i-1][1] + stepsizesim[i]*np.sin(thetasim[i])])
            counter+=1
        tracks.append(np.asarray(current_track))

    return tracks
