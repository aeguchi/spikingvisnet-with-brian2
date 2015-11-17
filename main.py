
from Parameters import *

#import imageImport as imimp
#import cv2

import netarch
import netplot
import GaborFilter

import numpy as np
import pylab as plt
import scipy

import glob
import os
import sys
import pickle
import errno

plotGabor = 0;
plotActivities = 0;
phases = [0,1,2] #testing only [0], testing and training [0,2]


try:
    os.makedirs("Results");
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise

try:
    os.makedirs("Results/"+experimentName);
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise

print "*** constructing the network ***"

vnet = netarch.visnet()
vplotter = netplot.plotter(vnet)
gf = GaborFilter.GaborFilter();

fileList_train = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_train.txt", dtype='str');
numObj_train = int(fileList_train[0]);
numTrans_train = int(fileList_train[1]);
fileList_test = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_test.txt", dtype='str');
numObj_test = int(fileList_test[0]);
numTrans_test = int(fileList_test[1]);

for phase in [0,1,2]:
    print "*** simulation phase: " + str(phase);
    if phase==1:
        vnet.setSynapticPlasticity(True)
        numObj = numObj_test;
        numTrans = numTrans_test;
    else:
        vnet.setSynapticPlasticity(False)
        numObj = numObj_train;
        numTrans = numTrans_train;
        
    # for img_fn in img_fns:
    FRrec = np.zeros((numObj_test,numTrans_test,nLayers+1, layerDim, layerDim));#1st layer is binding layer
    for index_obj in range(numObj):
        for index_trans in range(numTrans):
            print " ** Stimuli: " +str(phase) + " obj: " + str(index_obj) + ", trans: " + str(index_trans) +" **";
            index_img = index_obj * numTrans + index_trans;
            img_fn = os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/train/" + fileList_train[index_img + 2];
            # print img_fn
        
            # load gabor filtered ImageSurface
        
            #I = scipy.misc.imread(img_fn, flatten=True);
            inputImage = scipy.misc.imread(img_fn);
            if inputImage is None:
                print 'Failed to load image file:', img_fn
                sys.exit(1)

            #plt.imshow(I,cmap=plt.gray(),interpolation='none');
            #plt.colorbar();
            #plt.show()    

            resizedImage = gf.resizeImg(inputImage);
            
#             plt.imshow(I,cmap=plt.gray(),interpolation='none');
#             plt.colorbar();
#             plt.show()          
            
            res = gf.filtering(resizedImage);
            

            
            # filter the images and convert to normalised spikes
            
            res_norm = gf.scaleToUnit(res)
        
            vnet.setGaborFiringRates(res_norm)
        
            # run visnet simulation!
        
            print " ** running a simulation **"
            vnet.net.run(simulationTime * ms)
        
            # plot each image set
            if plotGabor:
                vplotter.plotGaborInput(inputImage, index_img, res, res_norm)
            FRrecTmp = np.zeros((nLayers+1, layerDim, layerDim));
            vplotter.plotLayers(inputImage, index_img, FRrecTmp)
            FRrec[index_obj,index_trans]=FRrecTmp;
            if (plotActivities):
                plt.show();
                
            if phase==0 or phase == 2:
                vnet.traceReset();
                #print "reset"

        vnet.traceReset();
    if phase==0:
        pickle.dump(FRrec, open("Results/"+experimentName+"/FR_0_blank.pkl", "wb"))
    elif phase==2:
        pickle.dump(FRrec, open("Results/"+experimentName+"/FR_1_trained.pkl", "wb"))
