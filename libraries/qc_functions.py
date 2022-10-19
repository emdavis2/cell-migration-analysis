import numpy as np
import matplotlib
import sys
import copy

def apply_qc(sample_in,minTrackLength,keep,trim,removemov,exclude):
    
    sample=copy.deepcopy(sample_in)
    #a sample containes movies, each movie contains tracks
    
    #initialize the status of all tracks in the sample as 1 (active)
    sampTrStatus=[]
    for imov in range(len(sample)):
        sampTrStatus.append([])
        for itr in range(len(sample[imov])):
            sampTrStatus[imov].append(1)
    
    
    #CONDITIONS ON TRACKS
    for imov in range(len(sample)):
        for itr in range(len(sample[imov])):        
            #if track length less than min length deactivate track
            if len(sample[imov][itr]) < minTrackLength:
                sampTrStatus[imov][itr]=0
    
    #exclude tracks after visual inspection
    #exclude=[[2,3],[2,4]]
    #exclude=[[1,7]]
                
    for i in range(len(exclude)):
        sampTrStatus[exclude[i][0]-1][exclude[i][1]-1]=0
    
    #only keep certain tracks
    #input: list of elements of the form [movie,track1,track2,...]
    #keep=[[2,3,4],[4,1]]
    #keep=[]
    for i in keep:
        #for all the tracks in movie i[0]-1 turn off all the tracks
        for itracks in range(len(sampTrStatus[ i[0]-1 ])):        
            sampTrStatus[i[0]-1][itracks]=0
        #turn on the desired tracks
        for itracks in i[1:]:
            sampTrStatus[i[0]-1][itracks-1]=1
            
    #remove movies
    #movie 3, 4 there are frames where the segmentation misses the cell
    for mov in removemov:
        for itr in range(len(sampTrStatus[mov-1])): 
            sampTrStatus[mov-1][itr]=0            
            
    #trim tracks
    #input: list (trim) of elements with the form [movie, track, begging frame, end frame]
    #trim=[[7,1,1,53]]
    for i in trim:    
        #print(i)
        #get 'frame' column
        framec=sample[i[0]-1][i[1]-1]['frame']
        #print (framec)
        #get index of desired first frame
        firstframe = framec.iloc[0]
        lastframe = framec.iloc[-1]
        if i[2] < firstframe or i[2] > lastframe:
            print('in movie',i[0],'track',i[1],'beggining of trimming',i[2],'is out of range',firstframe,lastframe)
            sys.exit()
        if i[3] < firstframe or i[3] > lastframe:
            print('in movie',i[0],'track',i[1],'end of trimming',i[3],'is out of range',firstframe,lastframe)
            sys.exit()

        if i[2] >= i[3] :
            print('in movie',i[0],'track',i[1],'end of trimming',i[3],' is smaller or equal than beggining of trimmig',i[2])
            sys.exit()
            
        ifirstframe = framec[framec==i[2]].index[0]
        #get index of desired last frame
        ilastframe = framec[framec==i[3]].index[0]
        #trim track
        sample[i[0]-1][i[1]-1]=sample[i[0]-1][i[1]-1].loc[ifirstframe:ilastframe+1]    
    
    return sampTrStatus, sample

