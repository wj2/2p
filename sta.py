"""
Create spike-triggered average of images to make a movie
"""

import numpy as np, scipy.io as sio, matplotlib.pyplot as plt
import scipy,time,os
from plotting import get_xm

def ndinterp(x,xp,yp):
    """n-dimensional interpolation
    Same interface as np.interp, except yp can be n-dimensional"""
    xp,yp = np.asarray(xp),np.asarray(yp)
    ind = scipy.interp(x,xp,np.arange(len(xp)))
    intpart = np.int32(np.trunc(ind))
    fracpart = ind - intpart
    shape = (-1,) + (1,)*(yp.ndim-1)
    return (1-fracpart).reshape(*shape)*yp[intpart]+ fracpart.reshape(*shape)*yp[intpart+1]

def test_ndinterp():
    assert (ndinterp([25,35],[20,30,40],[[10],[11],[12]]) == np.array([[10.5],[11.5]])).all()
        
    
def animate_frames(imgs,ts,save=False):
    """imgs is a sequence of 2d arrays. save them to files in current directory."""
    plt.ion()
    for t,img in zip(ts,imgs):
        plt.clf()
        plt.imshow(img)
        plt.title("t = %.2f seconds"%t)
        plt.draw()
        if save: plt.savefig("img%.2f.png"%t)
        else: time.sleep(.5)

        # mencoder "mf://*.png" -mf type=png:w=800:h=600:fps=5 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o psth.mpeg

if __name__ == "__main__":
    MAT = sio.loadmat("/home/joschu/Data/calcium/2010_05_24_0005.mat")
    TS = MAT["spktime"].flatten()
    TF = MAT["t_f"].flatten()
    X,M = get_xm()
    
    DTS = np.arange(0,3.5,.1)
    STAS = np.array([ndinterp(TS+t_delay,TF,X).mean(axis=0)-M
                    for t_delay in DTS])
