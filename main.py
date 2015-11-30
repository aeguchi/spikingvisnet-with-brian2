# date: 24/11/15
# author: Akihiro Eguchi
# description: a main code to run spiking network model.

from Parameters import *

import netarch
import netplot
import GaborFilter
import InfoAnalysis

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
phases = [0,1,1,1,1,1,2] #0:testing before training 1:training 2:testing after training
#phases = [0,1,2]
#phases = [1];

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
ia = InfoAnalysis.InfoAnalysis()

fileList_train = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_train.txt", dtype='str');
numObj_train = int(fileList_train[0]);
numTrans_train = int(fileList_train[1]);
fileList_test = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_test.txt", dtype='str');
numObj_test = int(fileList_test[0]);
numTrans_test = int(fileList_test[1]);

for phase in phases:
    print "*** simulation phase: " + str(phase) + " ***";
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
            print " ** obj: " + str(index_obj) + ", trans: " + str(index_trans) +" **";
            index_img = index_obj * numTrans + index_trans;
            img_fn = os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/train/" + fileList_train[index_img + 2];
            # print img_fn
        
            # load gabor filtered ImageSurface
        
            #I = scipy.misc.imread(img_fn, flatten=True);
            inputImage = scipy.misc.imread(img_fn);
            if inputImage is None:
                print 'Failed to load image file:', img_fn
                sys.exit(1)
            elif len(inputImage.shape)==3:
                inputImage = np.mean(inputImage[:,:,0:2],2);

            #plt.imshow(I,cmap=plt.gray(),interpolation='none');
            #plt.colorbar();
            #plt.show()    

            resizedImage = gf.resizeImg(inputImage);
            #print resizedImage
            
#             plt.imshow(I,cmap=plt.gray(),interpolation='none');
#             plt.colorbar();
#             plt.show()          
            
            res = gf.filtering(resizedImage);
            
            #res = np.power(res,2);
            

            
            # filter the images and convert to normalised spikes
            
            res_norm = gf.scaleToUnit(res)
        
            vnet.setGaborFiringRates(res_norm)
        
            # run visnet simulation!
            if phase in [0,2]:
                simulationTime = testingTime;
            else:
                simulationTime = trainingTime;
            
            vnet.net.run(simulationTime * ms)
            
        
            # plot each image set
            if plotGabor:
                vplotter.plotGaborInput(inputImage, index_img, res, res_norm,simulationTime)
            FRrecTmp = np.zeros((nLayers+1, layerDim, layerDim));
            vplotter.plotLayers(inputImage, index_img, FRrecTmp,simulationTime)
            
            FRrec[index_obj,index_trans]=FRrecTmp;
            if (plotActivities):
                plt.show();
                
            if phase==0 or phase == 2:
                vnet.traceReset();
        if phase==1:
            vnet.traceReset();
        
    if phase==0:
        pickle.dump(FRrec, open("Results/"+experimentName+"/FR_0_blank.pkl", "wb"))
    elif phase==2:
        pickle.dump(FRrec, open("Results/"+experimentName+"/FR_1_trained.pkl", "wb"))
    print "*** DONE ***"
    
plt.clf();
ia.singleCellInfoAnalysis();