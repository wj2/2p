"""
grab file names --> to read in data
pixels per line and lines per row --> to check res
microns per pixel --> for scale
relative and absolute times --> not sure
number of frames --> verify grabbed all of them
"""

import BeautifulSoup as bs,numpy as np,os
from tifffile import TIFFfile

def photonXMLGrabber(path):
    pieces = str.split(path,'/')
    directory = '/'
    for i in range(len(pieces)-1):
        directory += pieces[i]+'/'
    
    soup = bs.BeautifulStoneSoup(open(path))
    xRes,yRes,zRes = int(soup.find('key',key="linesPerFrame")['value']),int(soup.find('key',key="pixelsPerLine")['value']),len(soup('frame'))
    xMic,yMic = float(soup.find('key',key="micronsPerPixel_XAxis")['value']),float(soup.find('key',key="micronsPerPixel_YAxis")['value'])
    xPosW,yPosW,zPosW = soup.findAll('key',key="positionCurrent_XAxis"),soup.findAll('key',key="positionCurrent_YAxis"),soup.findAll('key',key="positionCurrent_ZAxis")
    xPos,yPos,zChange = [],[],[]
    for i in range(len(xPos)):
        xPos.append(float(xPosW[i]['value'])),yPos.append(float(yPosW[i]['value']))
        if i != 0:
            zChange.append(float(xPosW[i])-float(xPosW[i-1]))
    xAvg,yAvg,avgSliceWid = np.mean(xPos),np.mean(yPos),np.mean(zChange)
    xChange,yChange = xPos - xAvg,yPos - yAvg
    files = soup.findAll('file')
    names,channels = [],[]
    for fil in files:
        names.append(str(fil['filename']))
        channels.append(int(fil['channel']))
    print 'There are '+str(len(set(channels)))+' channels which would you like to analyse? (probably enter 1)'
    chan = raw_input('> ')
    chan = int(chan)
    shots = len(names)/len(set(channels))
    names = names[((chan-1)*shots):(chan*shots)]
    names.sort()
    return names,xMic,yMic,directory

def tiffSeriesReader(directory,names):
    prevD = os.getcwd()
    os.chdir(directory)
    tifs = []
    for i in range(len(names)):
        tifs.append(TIFFfile(names[i]).asarray())
    os.chdir(prevD)
    return tifs

#what to do about the multiple channels.
        
            


