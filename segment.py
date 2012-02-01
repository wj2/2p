"""
Segmentation of images

Current process
- Upsample image by a factor of 2
- Threshold image using local and global criterion
    local: intensity is above pth percentile in a surrounding block
    global: intensity is above pth percentile for whole image
- binary opening: separates blobs
- binary erosion: shrinks them to kill small blobs and yield few connected pixels in the middle of each blob
- label the blobs
- watershed segmentation using labeled blobs as markers (i.e., centers of basins)
- get rid of the blobs that don't fit the size criteria
"""

import numpy as np, scipy.ndimage as ndi, scipy.ndimage.morphology as snm, matplotlib.pyplot as plt, scipy as sp
from scipy import stats
from math import ceil

def local_maxima(X_nm):
    """Return binary array indicating local maxima of 2d array"""
    N,M = X_nm.shape
    return (
        np.r_[X_nm[:-1,:] >= X_nm[1:,:],np.zeros((1,M),bool)] & 
        np.r_[np.zeros((1,M),bool),X_nm[1:,:] >= X_nm[:-1,:]] & 
        np.c_[X_nm[:,:-1] >= X_nm[:,1:],np.zeros((N,1),bool)] &
        np.c_[np.zeros((N,1),bool),X_nm[:,1:] > X_nm[:,:-1]])

def test_local_maxima():
    X_nm = np.array([[0, 1, 1, 0],
                     [1, 2, 2, 1],
                     [0, 1, 1, 0]])
    print local_maxima(X_nm)


def region_borders(X_nm):
    """X_nm is a label array. Return a mask indicating borders"""
    N,M = X_nm.shape
    return (
        np.r_[X_nm[:-1,:] != X_nm[1:,:],np.zeros((1,M),bool)] |
        np.r_[np.zeros((1,M),bool),X_nm[1:,:] != X_nm[:-1,:]] | 
        np.c_[X_nm[:,:-1] != X_nm[:,1:],np.zeros((N,1),bool)] |
        np.c_[np.zeros((N,1),bool),X_nm[:,1:] != X_nm[:,:-1]])

def overlay_labels(labels):
    """Draw borders of labeled regions on the current figure
    0 is background. 
    Sometimes it looks like edges are missing. That's just because they're
    thin and don't get rendered"""
    labels_scaled = rescaled(labels,0,256)
    labels_colored = plt.get_cmap("jet")(labels_scaled)
    border_mask = region_borders(labels) & (labels > 0)
    labels_colored[~border_mask,:,3] = 0 # set alpha to zero
    return labels_colored

def overlay_full(labels):
    labels_scaled = rescaled(labels,0,256)
    labels_colored = plt.get_cmap("jet")(labels_scaled)
    border_mask = (labels > 0)
    labels_colored[~border_mask,:,3] = 0 # set alpha to zero
    return labels_colored
        
    
def watershed_segment(M,xM=None,yM=None):
    """Use watershed segmentation on an array. 
    Return an array where regions are given different integer labels"""

    if xM != None and yM != None:
        sel = np.ones((int(ceil(23.9*xM)),int(ceil(23.9*yM)))) # for opening
        sel2 = np.ones((int(ceil(127.2*xM)),int(ceil(127.2*yM)))) # for local thresholding
        sel3 = np.ones((int(ceil(11.9*xM)),int(ceil(11.9*yM)))) # for erosion
        ma,mi =(44245.21*xM*yM),(316.037*xM*yM) 
    else:
        selD = np.array([int(M.shape[0]*.012),int(M.shape[1]*.012)])
        selD = np.where(selD!=0,selD,1)
    
        sel2D = np.array([int(M.shape[0]*.12),int(M.shape[1]*.12)])
        sel2D = np.where(sel2D!=0,sel2D,1)

        sel3D = np.array([int(M.shape[0]*.01),int(M.shape[1]*.01)])
        sel3D = np.where(sel3D!=0,sel3D,1)


        sel = np.ones(selD) # for opening
        sel2 = np.ones(sel2D) # for local thresholding
        sel3 = np.ones(sel3D) # for erosion
        ma,mi = (M.shape[0]*M.shape[1]*.0075),(M.shape[0]*M.shape[1]*.0003)

    # get a few points in the center of each blob
    
    # threshold
    bw = ((M>=ndi.percentile_filter(M,80,footprint=sel2)))
    #& (M>=stats.scoreatpercentile(M.flatten(),80)))

    # open and erode
    blobs = snm.binary_opening(bw,structure=sel)
    blobs = snm.binary_erosion(blobs,structure=sel3,iterations=2)
    
    # label
    labels,_ = ndi.label(blobs)
    labels[labels > 0] += 1
    labels[0,0] = 1

    # rescale and cast to int16, then use watershed
    #M2 = rescaled(M,0,65000).astype(np.uint16)
    #newlabels = ndi.watershed_ift(M2,labels)
    newlabels = labels
    
    # get rid of groups unless they have the right number of pixels

    counts = np.bincount(newlabels.flatten())
    old2new = np.arange(len(counts)) 
    old2new[(counts < int(mi)) | (counts > int(ma))] = 0
    newlabels = old2new[newlabels]

    return newlabels

def watershed_segment_2(M,click_coords):
    """Use watershed segmentation on an array. 
    Return an array where regions are given different integer labels"""
    
    # todo: choose these structures based on aspect ratio of M and input parameters
    sel = np.ones((4,10)) # for opening
    sel2 = np.ones((15,75)) # for local thresholding
    sel3 = np.ones((2,5)) # for erosion
    # get a few points in the center of each blob
    
    # threshold
    #bw = ((M>=ndi.percentile_filter(M,80,footprint=sel2)) & (M>=scoreatpercentile(M.flatten(),60)))
    
    score = stats.percentileofscore(M.flatten(),M[int(click_coords[0][1]),int(click_coords[0][0])])
    bw = (M>=stats.scoreatpercentile(M.flatten(),score))

    # open and erode
    #bools = sp.zeros((M.shape[0],M.shape[1]),int)
    #bools[int(click_coords[0]),int(click_coords[1])] = 1
    #blobs = sp.where(bools == 1,True,False)
    blobs = snm.binary_opening(bw,structure=sel)
    blobs = snm.binary_dilation(blobs,iterations=3)
    blobs = snm.binary_erosion(blobs,structure=sel3)
    
    
    # label
    labels,_ = ndi.label(blobs)
    labels[labels > 0] += 1
    #labels[0,0] = 1

    # rescale and cast to int16, then use watershed
    M2 = rescaled(M,0,65000).astype(np.uint16)
    newlabels = ndi.watershed_ift(M2,labels)
    
    # get rid of groups unless they have the right number of pixels
    counts = np.bincount(newlabels.flatten())
    old2new = np.arange(len(counts))
    old2new[(counts < 100) | (counts > 600)] = 0
    newlabels = old2new[newlabels]
    
    return newlabels

def rescaled(M,newmin,newmax):
    """linearly rescale data in M
    returns a copy where values vary between newmin and newmax"""
    mmin,mmax = M.min(),M.max()
    M2 = M.copy()
    M2 -= mmin
    M2 *= (newmax-newmin) / (mmax-mmin)
    M2 += newmin
    return M2
    

if __name__ == "__main__":
    from plotting import get_xm
    import matplotlib.pyplot as plt
    import registration_alt as reg
    plt.clf()
    M = reg.correct_jitter("/home/wjj/2pProject/sample_data/turboreg013.npy")
    # upsample by a factor of 2
    M = ndi.interpolation.zoom(M[1],2)
    plt.figure(1)
    plt.imshow(M)
    LABELS = watershed_segment(M)
    plt.figure(1)
    overlay_labels(LABELS)
    plt.show()
