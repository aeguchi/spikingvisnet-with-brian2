# date: 24/11/15
# author: Akihiro Eguchi
# description: a main code to run spiking network model.

#from Parameters import *

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


#phases = [0,1,1,1,1,1,2] #0:testing before training 1:training 2:testing after training
phases = [0,1,1,1,1,1,2] #0:testing before training 1:training 2:testing after training
#phases = [0,1,2]
#phases = [1];

def loadParams(borrowed_globals):
    globals().update(borrowed_globals);

    
def runSimulation():
    try:
        os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/Results");
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    try:
        os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName);
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    print "*** constructing the network ***"
    
    vnet = netarch.visnet(globals())
    vplotter = netplot.plotter(vnet,globals())
    gf = GaborFilter.GaborFilter(globals());
    ia = InfoAnalysis.InfoAnalysis(globals())
    
    fileList_train = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_train.txt", dtype='str');
    numObj_train = int(fileList_train[0]);
    numTrans_train = int(fileList_train[1]);
    fileList_test = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_test.txt", dtype='str');
    numObj_test = int(fileList_test[0]);
    numTrans_test = int(fileList_test[1]);
    
    #plt.ion();
    count = 0;
    timeBegin = 0;
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
        SpikeRec = [];#this eventually becomes a list of SpikeRec[obj][trans][layer] #1st bind, 2nd 1st layer, 3rd 2nd layer..
        for index_obj in range(numObj):
            SpikeRec.append([])
            for index_trans in range(numTrans):
                SpikeRec[index_obj].append([]);
                print " ** obj: " + str(index_obj) + ", trans: " + str(index_trans) +" **";
                index_img = index_obj * numTrans + index_trans;
                img_fn = os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/train/" + fileList_train[index_img + 2];
            
                inputImage = scipy.misc.imread(img_fn);
                if inputImage is None:
                    print 'Failed to load image file:', img_fn
                    sys.exit(1)
                elif len(inputImage.shape)==3:
                    inputImage = np.mean(inputImage[:,:,0:2],2);

                resizedImage = gf.resizeImg(inputImage);
                res = gf.filtering(resizedImage);
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
                    vplotter.plotGaborInput(inputImage, index_img, res, res_norm, timeBegin,simulationTime)
                if plotActivities:
                    vplotter.plotLayers(inputImage, index_img, timeBegin,simulationTime)
                
                #to later save FR to the file
                if phase==0 or phase==2:
                    FRrecTmp = np.zeros((nLayers+1, layerDim, layerDim));
                    FRrecTmp[0] = vnet.getFiringRateMap(layerDim,vnet.spkdetBindingLayer,timeBegin,simulationTime);
                    SpikeRec[index_obj][index_trans].append(vnet.spkdetBindingLayer.spike_trains());
                    for i in range(nLayers):
                        FRrecTmp[i+1]=vnet.getFiringRateMap(layerDim,vnet.spkdetLayers[i],timeBegin,simulationTime);                
                        SpikeRec[index_obj][index_trans].append(vnet.spkdetLayers[i].spike_trains());
                    FRrec[index_obj,index_trans]=FRrecTmp;
                
                
                if (plotActivities or plotGabor):
                    plt.show(0);
                    vplotter.saveFigs(plotActivities,plotGabor,os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/"+str(count)+"_p"+str(phase)+"_o"+str(index_obj)+"_t"+str(index_trans));
                    #plt.show(block=False);
                    #plt.draw();
                plt.clf();
                
                if phase==0 or phase == 2:
                    vnet.traceReset();
                    
                timeBegin += simulationTime;
            if phase==1:
                vnet.traceReset();
        if phase==0:
            pickle.dump(FRrec, open(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/FR_0_blank.pkl", "wb"))
            pickle.dump(SpikeRec, open(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/Spikes_0_blank.pkl", "wb"))
        elif phase==2:
            pickle.dump(FRrec, open(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/FR_1_trained.pkl", "wb"))
            pickle.dump(SpikeRec, open(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/Spikes_1_trained.pkl", "wb"))
        print "*** DONE ***"
        count+=1;
        
    plt.clf();
    ia.singleCellInfoAnalysis();