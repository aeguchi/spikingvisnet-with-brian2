
from Parameters import *

import imageImport as imimp
import cv2

import netarch
import netplot

import numpy as np
import pylab as plt

import glob
import os
import sys
import pickle
import errno

plotGabor = 0;
#plotLayer = 0;  # 0:each images, 1: at end
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

vnet = netarch.visnet()
vplotter = netplot.plotter(vnet)

### Training ###

# img_fns = sorted(glob.iglob(os.path.split(os.path.realpath(__file__))[0] + "/images/training/*.png"))
# 
# if len(img_fns)!= nStim*nTrans:
#     print 'Error: the number of images files does not match',len(img_fns);
#     sys.exit(1)




fileList_train = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_train.txt", dtype='str');
numObj_train = int(fileList_train[0]);
numTrans_train = int(fileList_train[1]);
fileList_test = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_test.txt", dtype='str');numObj_train = int(fileList_train[0]);
numObj_test = int(fileList_test[0]);
numTrans_test = int(fileList_test[1]);

for phase in [0,1,2]:
    if phase==1:
        vnet.setSynapticPlasticity(True)
    else:
        vnet.setSynapticPlasticity(False)
        
    # for img_fn in img_fns:
    FRrec = np.zeros((numObj_test,numTrans_test,nLayers+1, layerDim, layerDim));#1st layer is binding layer
    for index_obj in range(numObj_train):
        for index_trans in range(numTrans_train):
            print "phase: " +str(phase) + ", obj: " + str(index_obj) + ", trans: " + str(index_trans);
            index_img = index_obj * numTrans_train + index_trans;
            img_fn = os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/train/" + fileList_train[index_obj * numTrans_train + index_trans + 2];
            # print img_fn
        
            # load gabor filtered ImageSurface
        
            img = cv2.imread(img_fn)
        
            if img is None:
                print 'Failed to load image file:', img_fn
                sys.exit(1)
        
            # filter the images and convert to normalised spikes
            
            filters = imimp.build_filters()
            res = imimp.process(img, filters)
            res_norm = imimp.to_spikes(res)
        
            vnet.setGaborFiringRates(res_norm)
        
            # run visnet simulation!
        
            vnet.net.run(simulationTime * ms)
        
            # plot each image set
            if plotGabor:
                vplotter.plotGaborInput(img, index_img, res, res_norm)
            FRrecTmp = np.zeros((nLayers+1, layerDim, layerDim));
            vplotter.plotLayers(img, index_img, FRrecTmp)
            FRrec[index_obj,index_trans]=FRrecTmp;
            if (plotActivities):
                plt.show();
                
            if phase==0 or phase == 2:
                vnet.traceReset();
        
    
            # index_img += 1
        vnet.traceReset();
    if phase==0:
        pickle.dump(FRrec, open("Results/"+experimentName+"/FR_0_blank.pkl", "wb"))
    elif phase==2:
        pickle.dump(FRrec, open("Results/"+experimentName+"/FR_1_trained.pkl", "wb"))
