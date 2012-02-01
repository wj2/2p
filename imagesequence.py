"""
Load image sequence from tiff file
"""

import numpy as np
import os
from tifffile import TIFFfile


def tiff2npy(tiff_name,npy_name=None):
    """Load tiff file, and write it to a npy file
    if npy_name is not specified, just use base name of given file (blah.tif -> blah.npy)"""
    tf = TIFFfile(tiff_name)
    npy_name = npy_name or os.path.splitext(tiff_name)[0]
    np.save(npy_name,tf.asarray())
     
