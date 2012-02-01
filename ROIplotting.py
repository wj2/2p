import numpy as np, scipy as sp, matplotlib.pyplot as plt
from scipy import ndimage as ndi
from scipy.signal import *

def find_avg_at_t(data_at_t,ROI_bitmask):
    points = np.where(ROI_bitmask > 0)
    if ROI_bitmask.shape != data_at_t.shape:
        print 'dif size'
        x = float(ROI_bitmask.shape[0])/float(data_at_t.shape[0])
        points = points[0]/x,points[1]/x
    summ = 0
    for i in range(len(points[1])):
        summ += data_at_t[int(points[0][i]),int(points[1][i])]
    return summ/float(len(points[1]))

def find_avg_over_tc(data,ROI_bitmask):
    over_tc = []
    #if ROI_bitmask.shape != data[0].shape:
    #    ROI_bitmask = ndi.interpolation.zoom(ROI_bitmask,(float(data.shape[1])/float(ROI_bitmask.shape[0])))
    for t in range(data.shape[0]):
        over_tc.append(find_avg_at_t(data[t],ROI_bitmask))
    return over_tc

#def hpFilter(data,bitmask):

def lpFilter(data,bitmask,order,cFreq):
    fil = iirfilter(order,cFreq,btype='lowpass')
    return lfilter(fil[0],fil[1],find_avg_over_tc(data,bitmask)),fil
    
def spikeDet(data,bitmask,mult=.3,filterOrder=3,filterCritFreq=.05):
    train,fil = lpFilter(data,bitmask,filterOrder,filterCritFreq)
    #th = np.median(train)+mult*np.std(train)
    th = np.median(train)+mult*np.median(abs(train-np.median(train)))
    postrig = train > th
    return train,np.diff(postrig).nonzero()[0][::2],fil            

def spikeDetL(data,bitmask,fil,mult=.3):
    train = lfilter(fil[0],fil[1],find_avg_over_tc(data,bitmask))
    #th = np.median(train)+mult*np.std(train)
    th = np.median(train)+mult*np.median(abs(train-np.median(train)))
    postrig = train > th
    return train,np.diff(postrig).nonzero()[0][::2]

def plot_it(data,ROI_bitmask,noise_reduc=None):
    data_over_tc = find_avg_over_tc(data,ROI_bitmask)
    if noise_reduc != None:
        data_over_tc = noise_reduction(data_over_tc,noise_reduc)
    plt.plot(data_over_tc)
    return None
    
#here lie shenanigans
def lpFilterS(fluor,order,cFreq):
    fil = iirfilter(order,cFreq,btype='lowpass')
    return lfilter(fil[0],fil[1],fluor),fil

def spikeDetS(fluor,mult=.3,filterOrder=3,filterCritFreq=.05):
    train,fil = lpFilterS(fluor,filterOrder,filterCritFreq)
    #th = np.median(train)+mult*np.std(train)
    th = np.median(train)+mult*np.median(abs(train-np.median(train)))
    postrig = train > th
    return train
    
    
