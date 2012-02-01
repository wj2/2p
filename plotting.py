import matplotlib.pyplot as plt, numpy as np, scipy.ndimage as ndi, scipy as sp,sys, imagesequence as ims
from matplotlib.widgets import SpanSelector,Button
import registration as reg, segment as seg,ROIplotting as roip
from tifffile import TIFFfile
import xmlParser as xml
from scipy.stats import pearsonr


def get_xm(path):
    """load data and calculate the mean image"""
    X = np.load(path)
    M = np.array(np.mean(X,axis=0))
    return X,M   


def find_points(watershed_label):
    """finds the points in all individual ROIs of the watershed algorithm output"""
    ROIs = []
    for i in xrange(1,watershed_label.max()+1):
        instances = np.where(watershed_label == i)
        ROI = []
        for j in xrange(len(instances[1])):
            ROI.append([instances[0][j],instances[1][j]])
        if len(ROI)>1:
            ROIs.append(ROI)
    return ROIs


def create_bitmask(ROI,M):
    """creates mask where 1 values are points in the region and 0 values are points outside the region"""
    mask = np.zeros((M.shape[0],M.shape[1]))
    for roi in ROI:
        #print ROI[i][0],ROI[i][1]
        mask[roi[0],roi[1]] = 1
    return mask


def match_point(click,bitmasks):
    """matches point from user input to a specific bitmask"""
    for mask in bitmasks:
        if mask[int(click[0][1]),int(click[0][0])] > 0:
            return i


def delete_roi_2D(click,bitmasks):
    """deletes an ROI based on user input"""
    for i,mask in enumerate(bitmasks):
        if mask[int(click[0][1]),int(click[0][0])] > 0:
            return np.delete(bitmasks,i,0)

def loadData(path):
    spl = str.split(path,'.')
    if spl[-1] == 'npy':
        data = np.load(path)
    if spl[-1] == 'tif':
        data = TIFFfile(path).asarray()
    if spl[-1] == 'xml':
        names,xMix,yMic,dire = xml.photonXMLGrabber(path)
        data = xml.tiffSeriesReader(dire,names)
    return data


#def masks2Labels(masks):
    """converts masks back to labels, not currently used"""
    """
    overall = np.zeros((masks.shape[0],masks.shape[1]))
    for i in range(len(masks)):
        points = np.where(masks[i] == 0)
        for point in points:
            overall[point[0],point[1]] = i+1
    return overall
    """

def deleteMultiDetection(bitmasks,data):
    """finds unique ROIs in three dimensions, returns a list of indices of the principal z-slice of
    all ROIs as well as changes the bitmasks such that points not in the region are zero and points
    in the region are the same non-zero as all points in the ROI across the z-axis.

    for example, if a and b are ROI on different but adjacent z-slices and their centres are close
    enough such that the algorithm believes them to be duplicates of the same cell body, both arrays
    will have the same integer value in all points considered in the region."""
    prinROI,saveROI = [],[]
    count = 1
    for i in xrange(len(bitmasks)):
        for j in xrange(len(bitmasks[i])):
            if len(np.where(bitmasks[i][j] == 1)[0]) != 0:
                a = np.where(bitmasks[i][j] == 1)
                xavg1,yavg1 = np.mean(a[0]),np.mean(a[1])
                sameROILoc = [[i,j]]
                avg = [roip.find_avg_at_t(data[i],bitmasks[i][j])]
                for k in xrange(i+1,len(bitmasks)-1):
                    q = 0
                    for l in xrange(len(bitmasks[k])):
                        b = np.where(bitmasks[k][l] >= 1)
                        xavg2,yavg2 = np.mean(b[0]),np.mean(b[1])
                        if abs(xavg1 - xavg2) < bitmasks[i][j].shape[0]*.1 and abs(yavg1-yavg2) < bitmasks[i][j].shape[1]*.1:
                            sameROILoc.append([k,l])
                            avg.append(roip.find_avg_at_t(data[k],bitmasks[k][l]))
                            q = 1
                    if q == 0:
                    #saveROI.append(sameROILoc[np.argmax(avg)])
                    #sameROILoc = np.delete(sameROILoc,np.argmax(avg),0)
                        count += 1
                        prinROI.append(sameROILoc[np.argmax(avg)])
                        for z in xrange(len(sameROILoc)):
                            bitmasks[sameROILoc[z][0]][sameROILoc[z][1]] = np.where(bitmasks[sameROILoc[z][0]][sameROILoc[z][1]] >= 1,count,0)
                        break
    prinROI = [list(ro) for ro in set(tuple(roi) for roi in prinROI)]
    prinROI.sort()
    return bitmasks,prinROI



def gtCompare(gtSpikes,data,bitmasks):
    avg,avgB,u = [],[],0
    for spk in gtSpikes:
        for i,mask in enumerate(bitmasks):
            if u == 0:
                avg.append(0.0)
                avgB.append(np.median(roip.find_avg_over_tc(data,mask)))
            avg[i] += (sum(roip.find_avg_over_tc(data[spk:spk+10],mask))/10.0)-avgB[i]
        u = 1
    return avg

def corCompare(fluor,data,bitmasks):
    r = []
    for mask in bitmasks:
        tc = roip.find_avg_over_tc(data,mask)
        r.append(pearsonr(fluor,tc))
    return r
        
            

#def deleteMultiDetection_backup(bitmasks,data):
    """back up of a different way of dealing with multiple detections of the same physical ROI across z."""
    """
    princRIO,saveROI = [],[]
    for i in range(len(bitmasks)):
        for j in range(len(bitmasks[i])):
            a = np.where(bitmasks[i][j] == 0)
            xavg1,yavg1 = np.mean(a[0]),np.mean(a[1])
            sameROILoc = [[i,j]]
            avg = [roip.find_avg_at_t(data[i],bitmasks[i][j])]
            for k in range(i+1,len(bitmasks)-1):
                q = 0
                for l in range(len(bitmasks[k])):
                    b = np.where(bitmasks[k][l] == 0)
                    xavg2,yavg2 = np.mean(b[0]),np.mean(b[1])
                    if abs(xavg1 - xavg2) < bitmasks[i][j].shape[0]*.1 and abs(yavg1-yavg2) < bitmasks[i][j].shape[1]*.1:
                        sameROILoc.append([k,l])
                        avg.append(roip.find_avg_at_t(data[k],bitmasks[k][l]))
                        q = 1
                if q == 0:
                    #saveROI.append(sameROILoc[np.argmax(avg)])
                    princROI.append(sameROILoc[np.argmax(avg)])
                    sameROILoc = np.delete(sameROILoc,np.argmax(avg),0)
                    for z in range(len(sameROILoc)):
                        saveROI.append(sameROILoc[z])
                    q = 1
                    break

    saveROI = [list(ro) for ro in set(tuple(roi) for roi in saveROI)]
    saveROI.sort()
    saveROI.reverse()
    #return saveROI
    for roi in saveROI:
        bitmasks[roi[0]] = np.delete(bitmasks[roi[0]],roi[1],0)
    return bitmasks
    """
"""    
def process2D(path):
    (X,M) = get_xm(path)
    M = ndi.interpolation.zoom(M,2)
    plt.figure(1)
    plt.imshow(M,aspect='equal',cmap=plt.cm.gray)
    LABELS = seg.watershed_segment(M)
    seg.overlay_labels(LABELS)
    print 'Finding the points in each detected ROI...'
    ROIs = find_points(LABELS)
    print 'Done.'
    print 'Would you like to select more ROIs?'
    cont = raw_input('> ')
    while cont == 'y':
        print 'Click in the upper left and then the lower right of a region that has an undetected neuron.'
        click_coords = plt.ginput(1,timeout=300)
        print click_coords
        #self.more_labels = local_watershed(self.click_coords,self.M)
        #self.more_labels = detect_roi_inrgn(self.click_coords,self.M)
        more_labels = seg.watershed_segment_2(M,click_coords)
        plt.figure(1)
        seg.overlay_labels(more_labels)
        x = find_points(more_labels)
        print 'We found '+str(len(x))+' new ROIs.'
        print 'Did that work? (y/n)'
        work = raw_input('> ')
        if work == 'y':
            ROIs.append(x)
        print 'Would you like to select more ROIs? (y/n)'
        cont = raw_input('> ')
    masks = np.ones((len(ROIs),M.shape[0],M.shape[1]))
    for k in range(len(ROIs)):
        masks[k] = create_bitmask(ROIs[k],M)
    print 'Save?'
    sav = raw_input('> ')
    if sav == 'y':
        np.save('ROI_bitmasks_for_'+name,masks)
    print 'Would you like to look at data from ROIs?'
    cont2 = raw_input('> ')
    fig = 1
    while cont2 == 'y':
        click_coords = plt.ginput(1,timeout=300)
        mask = match_point(click_coords,masks)
        fig += 1
        plt.figure(fig)
        if mask != None:
            roip.plot_it(X,mask,noise_reduc=20)
        else:
            print 'You missed.'   
    print 'More? (y/n)'
    cont2 = raw_input('> ')
    return None
"""
        
def process3D(data,xM=None,yM=None):
    """used to apply watershed to 3D or even 2D data"""
    #print 'Would you like to upsample by a factor of two? This will take longer, but might improve the result. (y/n)'
    #yn = raw_input('> ')
    allmasks,alllabels = [],[]
    for i in range(len(data)):
        #if yn == 'y':
        #    data[i] = ndi.interpolation.zoom(data[i],2)
        alllabels.append(seg.watershed_segment(data[i],xM,yM))
        ROIs = find_points(alllabels[i])
        data = np.array(data)
        masks = np.ones((len(ROIs),data.shape[1],data.shape[2]))
        for j in range(len(ROIs)):
            masks[j] = create_bitmask(ROIs[j],data[i])
        allmasks.append(masks)
    allmasks,prinROI = deleteMultiDetection(allmasks,data)
    return allmasks,prinROI

def process2D(data,xM=None,yM=None):
    if len(data.shape) > 2:
       _ ,data = reg.correct_jitter(data)
    alllabels = seg.watershed_segment(data,xM,yM)
    ROIs = find_points(alllabels)
    masks = np.ones((len(ROIs),data.shape[0],data.shape[1]))
    for i,roi in enumerate(ROIs):
        masks[i] = create_bitmask(roi,data)
    return masks,[],data

"""
def findROIsFrom2Photon(path,dim='3D'):
    name = path.split('/')
    name = name[len(name)-1]
    name,_ = name.split('.')
    if dim == '3D':
        process3D(path)
    if dim == '2D':
        process2D(path)
    else:
        print dim+' not supported.'
    return None
"""


COLORS = ["r","y","g","b"]*4
class CaViewer(object):
    """
    Class for the GUI. Takes single-file input in .tif or .npy, can take a series of .tif files if the
    .xml metadata file is provided in the same folder as the .tif files. 
    
    This is uses a combination of set_data and imshow because set_data seems to completely reset the image, 
    and did not work for overlaying the ROI labels. 
    
    The ROI has five buttons, two for scrolling forward and backward, one for deleting ROIs currently being
    used, one for designating an ROI that should be used, and, finally, one to return the ROIs presently 
    being used into a npy file.
    """
    def __init__(self,x = None, m = None):
        xMic,yMic = None,None
        
        try:
            path = sys.argv[1]
        except:
            print 'Please enter the path to the .xml metadata of your data or the single .tif or .npy file.'
            path = raw_input('> ')
        #handles multiple input types
        self.data = loadData(path)

        
        #self.data = np.load('sample_data/3Dtiffdata.npy')
        #processes the data, this step takes longest
        try:
            typ = sys.argv[2]
        except:
            print 'Is the data 3D or 2D?'
            typ = raw_input('> ')
        if typ == None or typ == '3D' or typ == '3d':
            self.bitmasks,self.prinROI = process3D(self.data,xMic,yMic)
            typ = 3
        elif typ == '2D' or typ == '2d':
            self.bitmasks,self.prinROI,self.avg = process2D(self.data,xMic,yMic)
            self.tc_data = self.data
            typ = 2
        self.t = 0  
        self.fig = plt.figure(1)
        self.lenData,self.i = len(self.data),0
        

        def loadI(i,bitmasks,prinROI):
            #refreshes the display on the GUI, declared here due to use in initialisation
            if typ == 3:
            #self.im.set_data(self.data[i])
                self.ax_img.imshow(self.data[i],aspect='equal',cmap=plt.cm.gray,interpolation='nearest')
                for j in range(len(bitmasks[i])):
                    self.ax_img.imshow(seg.overlay_labels(bitmasks[i][j]),aspect='equal',interpolation='nearest')
                    for prin in prinROI:
                        if prin == [i,j]:
                            self.ax_img.imshow(seg.overlay_full(bitmasks[i][j]),aspect='equal',interpolation='nearest')
            if typ == 2:
                self.ax_img.imshow(self.avg,aspect='equal',cmap=plt.cm.gray,interpolation='nearest')
                for j in range(len(bitmasks)):
                    self.ax_img.imshow(seg.overlay_labels(bitmasks[j]),aspect='equal',interpolation='nearest')
            
            """
            for j in range(len(bitmasks[i])):
                self.im.set_data(seg.overlay_labels(bitmasks[i][j]))
                for prin in prinROI:
                    if prin == [i,j]:
                        self.im.set_data(seg.overlay_full(bitmasks[i][j]))
            """


        self.ax_img = self.fig.add_axes([.1,.40,.8,.55])
        self.ax_plot = self.fig.add_axes([.1,.1,.8,.25])
        self.ax_but1 = self.fig.add_axes([.1,.05,.1,.04])
        self.ax_but2 = self.fig.add_axes([.2,.05,.1,.04])
        self.ax_but3 = self.fig.add_axes([.3,.05,.1,.04])
        self.ax_but4 = self.fig.add_axes([.4,.05,.1,.04])
        self.ax_but5 = self.fig.add_axes([.8,.05,.1,.04])
        self.ax_but6 = self.fig.add_axes([.7,.05,.1,.04])
        self.ax_but7 = self.fig.add_axes([.6,.05,.1,.04])
        self.ax_but8 = self.fig.add_axes([.5,.05,.1,.04])
        self.ax_but9 = self.fig.add_axes([.8,.01,.1,.04])
        #self.ax_fullrange = self.fig.add_axes([.1,.275,.8,.15])
        #self.ax_zoom = self.fig.add_axes([.1,.1,.8,.15])
        if typ == 2:
            self.im = self.ax_img.imshow(self.avg,aspect='equal',cmap=plt.cm.gray,interpolation='nearest')
            for j in range(len(self.bitmasks)):
                self.ax_img.imshow(seg.overlay_labels(self.bitmasks[j]),aspect='equal',interpolation='nearest')
        if typ == 3:
            self.im = self.ax_img.imshow(self.data[self.i],aspect='equal',cmap=plt.cm.gray,interpolation='nearest')
            loadI(self.i,self.bitmasks,self.prinROI)
        #self.LABELS = seg.watershed_segment(self.M)
        #span = SpanSelector(self.ax_fullrange, self.onselect, 'horizontal', useblit=True,rectprops=dict(alpha=0.5, facecolor='red') )
        #self.line1, = self.ax_fullrange.plot(self.x,self.y)
        #self.line2, = self.ax_zoom.plot(self.x,self.y)
        self.but1 = Button(self.ax_but1,"<---")
        self.but2 = Button(self.ax_but2,"DELETE")
        self.but3 = Button(self.ax_but3,"USE")
        self.but4 = Button(self.ax_but4,"--->")
        self.but5 = Button(self.ax_but5,"SHOW SPIKES")
        self.but6 = Button(self.ax_but6,"CLEAR PLOT")
        self.but7 = Button(self.ax_but7,"ALL SPIKES")
        self.but8 = Button(self.ax_but8,"SINGLE SPIKE")
        self.but9 = Button(self.ax_but9,"SAVE MASKS")

        def onscroll(self,event):
            if event.button == 'up':
                prevS(event)
            else:
                nextS(event)
        
        #self.ax_img.set_yticks([])
        #self.ax_img.set_xticks([])
        #self.but3 = Button(self.ax_but3,
        self.roi_pts = None
        
        def prevS(event):
            if self.i-1 >= 0:
                self.i-=1
                print str(self.i)+'/'+str(self.lenData)
                loadI(self.i,self.bitmasks,self.prinROI)

        def nextS(event):
            if self.i+1 <= self.lenData:
                self.i+=1
                print str(self.i)+'/'+str(self.lenData)
                loadI(self.i,self.bitmasks,self.prinROI)

        def delROI(event):
            click = plt.ginput(1,300)
            if typ == 2:
                self.bitmasks = delete_roi_2D(click,self.bitmasks)
            if typ == 3:
                j = match_point(click,self.bitmasks)
                self.prinROI.remove([self.i,j])
            loadI(self.i,self.bitmasks,self.prinROI)
            print 'ROI will not be used.'

        def useROI(event):
            click = plt.ginput(1,300)
            j = match_point(click,self.bitmasks[self.i])
            self.prinROI.append([self.i,j])
            self.prinROI.sort()
            loadI(self.i,self.bitmasks,self.prinROI)
            #self.uniqueMasks[i] = vstack([self.bitmasks[si],addROI[newaxis,...]])
            print 'ROI will be used.'
            
        def returnROIs(event):
            masks = []
            for i in range(len(self.bitmasks)):
                for roi in self.prinROI:
                    if roi[0] == i:
                        masks.append(bitmasks[i][roi[1]])
            print 'What would you like to name this dataset?'
            name = raw_input('> ')
            np.save(name,masks)
            print 'ROIs have been saved.'

        def showSpikes(event):
            click = plt.ginput(1,300)
            if typ == 2:
                j = match_point(click,self.bitmasks)
                if self.t == 0:
                    train,spikes,self.fil = roip.spikeDet(self.tc_data,self.bitmasks[j])
                    self.t = 1
                elif self.t == 1:
                    train,spikes = roip.spikeDetL(self.tc_data,self.bitmasks[j],self.fil)
                else:
                    print 'Something is wrong.'
            #self.plo = plt.figure(2)
                self.ax_plot.plot(train)
                self.ax_plot.plot(spikes,np.ones(np.size(spikes))*35000,'r.')
                print 'Spiketimes for ROI '+str(j)+': ',spikes
            if typ == 3:
                print 'Not currently supported.'

        def clearPlot(event):
            self.ax_plot.clear()

        def saveAllSpkTrains(event):
            arch = []
            if typ == 2:
                for j in range(len(self.bitmasks)):
                    print 'Mask '+str(j+1)+' of '+str(len(self.bitmasks))
                    if self.t == 0:
                        train,spikes,self.fil = roip.spikeDet(self.tc_data,self.bitmasks[j])
                        self.t = 1
                    else:
                        train,spikes = roip.spikeDetL(self.tc_data,self.bitmasks[j],self.fil)
                    arch.append([train,spikes])
            elif typ == 3:
                for roi in self.prinROI:
                    if self.t == 0:
                        train,spikes,self.fil = roip.spikeDet(self.tc_data,self.bitmasks[roi[0]][roi[1]])
                    else:
                        train,spikes = roip.spikeDetL(self.tc_data,bitmasks[roi[0]][roi[1]],self.fil)
                        self.t = 1
                    arch.append([roi,train,spikes])
            print 'What would you like to name this dataset?'
            name = raw_input('> ')
            np.save(name,arch)

        def saveSingleSpkTrain(event):
            click = plt.ginput(1,300)
            if typ == 2:
                j = match_point(click,self.bitmasks)
                if self.t == 0:
                    train,spikes,self.fil = roip.spikeDet(self.tc_data,self.bitmasks[j])
                    self.t = 1
                elif self.t == 1:
                    train,spikes = roip.spikeDetL(self.tc_data,self.bitmasks[j],self.fil)
                arch = [train,spikes]
                print 'What would you like to name the data from this ROI?'
                name = raw_input('> ')
                np.save(name,arch)
                                                 
                    
                
            

        self.but1.on_clicked(prevS)
        self.but2.on_clicked(delROI)
        self.but3.on_clicked(useROI)
        self.but4.on_clicked(nextS)
        self.but5.on_clicked(showSpikes)
        self.but6.on_clicked(clearPlot)
        self.but7.on_clicked(saveAllSpkTrains)
        self.but8.on_clicked(saveSingleSpkTrain)
        self.but9.on_clicked(returnROIs)

        self.fig.canvas.draw()
    
        
if __name__ == "__main__":
    #fig.canvas.mpl_connext('scroll_event',c.onscroll)
    c = CaViewer()
    plt.show()
    pass
