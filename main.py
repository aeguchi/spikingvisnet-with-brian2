
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

plotGabor = 1;
plotLayer = 0; # 0:each images, 1: at end

vnet = netarch.visnet()
vplotter = netplot.plotter(vnet)

### Training ###

# img_fns = sorted(glob.iglob(os.path.split(os.path.realpath(__file__))[0] + "/images/training/*.png"))
# 
# if len(img_fns)!= nStim*nTrans:
#     print 'Error: the number of images files does not match',len(img_fns);
#     sys.exit(1)

fileList_train = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_train.txt",dtype='str');
numObj_train = int(fileList_train[0]);
numTrans_train = int(fileList_train[1]);
fileList_test = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_test.txt",dtype='str');numObj_train = int(fileList_train[0]);
numTrans_test = int(fileList_test[0]);
numTrans_test = int(fileList_test[1]);


#index_img=0

vnet.setSynapticPlasticity(True)

#for img_fn in img_fns:
for index_obj in range(numObj_train):
    for index_trans in range(numTrans_train):
        index_img = index_obj*numTrans_train + index_trans;
        img_fn = os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/train/" + fileList_train[index_obj*numTrans_train + index_trans + 2];
        #print img_fn
    
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
    
        vnet.net.run(simulationTime*ms)
    
        # plot each image set
    
        vplotter.plotGaborInput(img, index_img, res, res_norm)
        vplotter.plotLayers(img, index_img, plotLayer=plotLayer)
    
        #index_img += 1
    
vnet.setSynapticPlasticity(False)

#for img_fn in img_fns:
for index_obj in range(numObj_test):
    for index_trans in range(numTrans_test):
        index_img = index_obj*numTrans_test + index_trans;
        img_fn = os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/test/" + fileList_test[index_obj*numTrans_train + index_trans + 2];
        #print img_fn
    
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
    
        vnet.net.run(simulationTime*ms)
    
        # plot each image set
    
        if plotGabor:
            vplotter.plotGaborInput(img, index_img, res, res_norm)
        vplotter.plotLayers(img, index_img, plotLayer=plotLayer)
    
        #index_img += 1
