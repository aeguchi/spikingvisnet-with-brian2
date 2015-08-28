
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
plotLayer = 0; #0:each images, 1: at end

vnet = netarch.visnet()
vplotter = netplot.plotter(vnet)

### Training ###

img_fns = sorted(glob.iglob(os.path.split(os.path.realpath(__file__))[0] + "/images/training/*.png"))

if len(img_fns)!= nStim*nTrans:
    print 'Error: the number of images files does not match',len(img_fns);
    sys.exit(1)

index_img=0;

for img_fn in img_fns:

    #print img_fn
    #load gabor filtered ImageSurface

    img = cv2.imread(img_fn)

    if img is None:
        print 'Failed to load image file:', img_fn
        sys.exit(1)
    
    filters = imimp.build_filters()
    res = imimp.process(img, filters)
    res_norm = imimp.to_spikes(res)

    for index_filter in range(0,len(thetaList)):
        r = np.reshape(np.mean(res_norm[index_filter],axis=2),(layerGDim*layerGDim));
        #print r
        vnet.layerG[index_filter].rates= r * Rmax;    #To be fixed
        #print vnet.layerG[index_filter].rates

    # run visnet simulation!

    vnet.net.run(simulationTime*ms)

    vplotter.plotGaborInput(img, res, res_norm)
    vplotter.plotLayers(img, plotLayer=0)

    index_img += 1
