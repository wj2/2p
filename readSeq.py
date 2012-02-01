import os
import numpy as np, imagesequence as ims

def readTiffSeq(direc,baseName):
    x = os.getcwd()
    os.chdir(direc)
    files = filter(lambda x: baseName in x and '.tif' in x,os.listdir(os.getcwd()))
    files = np.sort(files)
    for i in range(len(files)):
        ims.tiff2npy(files[i])
    newfiles = filter(lambda x: baseName in x and '.npy' in x,os.listdir(os.getcwd()))
    newfiles = np.sort(newfiles)
    shapedata = np.load(newfiles[0])
    overall = np.ones((len(newfiles),shapedata.shape[0],shapedata.shape[1]))
    for j in range(len(newfiles)):
        overall[j] = np.load(newfiles[j])
    np.save(baseName,overall)
    os.chdir(x)
    return overall
 
