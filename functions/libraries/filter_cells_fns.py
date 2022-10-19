import numpy as np


def remove_multiple_nuclei_cells(labeledcellsmask,labelednucsmask):
    out = np.copy(labeledcellsmask)
    labels=set(out.ravel())
    for i in labels:
        if i!=0: #ignore background           
            #get obect i coordinates 
            icoord=np.where(out==i)
            #get values of nuclear mask for those coordinates
            nucvalues=labelednucsmask[icoord]
            #get labels list (get individual elements)
            nuclabels=set(nucvalues)
            #if there is more than one nuclei (2 objects counting backround)
            if len(nuclabels) > 2: 
                #remove cell object:
                out[out==i]=0
    return out
                
        
def remove_large_objects(labeledmask,maxarea):
    out = np.copy(labeledmask)
    component_sizes = np.bincount(labeledmask.ravel())
    too_large = component_sizes > maxarea
    too_large_mask = too_large[labeledmask]
    out[too_large_mask] = 0
    return out

def remove_touching_edge(labeledmask):    
    out = np.copy(labeledmask)
    #scan each edge of the image
    width=out.shape[1]
    height=out.shape[0]    
    j=0
    #zeroloc=np.where(out[:,j] != 0)
    #for k in zeroloc:
    #    out[out==out[zeroloc,j]]=0
        
    for i in range(height):
        if out[i,j] != 0:
            out[out==out[i,j]]=0
    j=width-1
    for i in range(height):
        if out[i,j] != 0:
            out[out==out[i,j]]=0
    i=0
    for j in range(width):
        if out[i,j] != 0:
            out[out==out[i,j]]=0
    i=height-1
    for j in range(width):
        if out[i,j] != 0:
            out[out==out[i,j]]=0    
    return out