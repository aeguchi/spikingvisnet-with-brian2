
from imageImport import *
import simplejson as json
import serial_json2
from PlotSim import *
import sys
import  glob
#import nest
import os

#from nest import raster_plot
#import nest.topology as topp




#http://brian2.readthedocs.org/en/latest/resources/tutorials/1-intro-to-brian-neurons.html

def RunSim():
    
    index_img=0;
    for img_fn in img_fns:
        #print img_fn
        #load gabor filtered ImageSurface
        img = cv2.imread(img_fn)
        if img is None:
            print 'Failed to load image file:', img_fn
            sys.exit(1)
        
        filters = build_filters()
        res = process(img, filters)
        
        #convert from gabor filtered inputs to spikes
        #normalize
        res_norm=res/numpy.max(res);
        res_norm=1-res_norm;
        
        for index_filter in range(0,len(thetaList)):
            r = numpy.reshape(numpy.mean(res_norm[index_filter],axis=2),(layerGDim*layerGDim));
            #print r
            layerG[index_filter].rates= r * Rmax;    #To be fixed
        #print layerG[index_filter].rates
        
        net.run(simulationTime*ms)
        
        store('trained')
        
        
        if (plotGabor == 1 ):
            PlotGabor(img,res_norm,res,index_img)
        if (plotLayer == 1 ):
            PlotLayer(img,res_norm,res,index_img, index_filter)
        
        index_img+=1;
    
    
    
    return 0


### Training ###
#trainingImages = sorted(glob.iglob("/images/training/*.png"))
trainingImages = sorted(glob.iglob(os.path.split(os.path.realpath(__file__))[0] + "/images/test/*.png"))

store('initialized')    #randomised network

if not os.path.exists(os.path.split(os.path.realpath(__file__))[0] + "/randomised/"):
    os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/randomised/")
if not os.path.exists(os.path.split(os.path.realpath(__file__))[0] + "/trained/"):
    os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/trained/")
if not os.path.exists(os.path.split(os.path.realpath(__file__))[0] + "/FR/"):
    os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/FR/")


img_fns = []
for img_fn in trainingImages:
    img_fns.append(img_fn)

if len(img_fns)!=nStim*nTrans:
    print 'Error: the number of images files does not match',len(img_fns);
    sys.exit(1)

for trial in range(trialNb) :

    print('Trial '+ str(trial+1) + ' out of '+str(trialNb) )
    restore('initialized')
    
    plasticON(0)
    RunSim()
    
    for layer in range(0,nLayers):
        filep= os.path.split(os.path.realpath(__file__))[0] + '/FR/' + 'FR_b_Ex'+str(layer)+'_trial'+str(trial)
        numpy.save(filep,spkdetLayers[layer].spike_trains())
        #with open(filep,'w') as FR:
            #json.dump( spkdetLayers[layer].spike_trains(), FR)
    
    
    plasticON(1)
    RunSim()
    
    for test in range(testNb):
        print('Test '+ str(test+1) + ' out of '+str(testNb) )
        restore('trained')
        plasticON(0)
        RunSim()
        
       
        for layer in range(0,nLayers):
            filep= os.path.split(os.path.realpath(__file__))[0] + '/FR/' + 'FR_t_Ex'+str(layer)+'_trial'+str(trial)+'_test'+str(test)
            numpy.save(filep,spkdetLayers[layer].spike_trains())
#            with open(filep,'w') as FR:
#               json.dump(spkdetLayers[layer].spike_trains(), FR)


    if (plotLayer == 2 ):
        PlotLayer(img,res_norm,res,index_img,index_filter)

#To-DO: trainNetworkWith
    
#Testing



