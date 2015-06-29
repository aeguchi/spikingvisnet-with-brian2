from Parameters import *
from imageImport import *
import sys
import numpy,pylab as plt, glob
import nest
from random import random



### Constructing Network ###
nest.Models()
ndict = {"I_e": 180.0, "tau_m": 20.0}
nest.SetDefaults("iaf_neuron", ndict)
inputNeurons = nest.Create("iaf_neuron", 10*10);
Vth=-55.                  
Vrest=-70.  
spikedetectors = nest.Create("spike_detector", 10*10, params={"withgid": True, "withtime": True})
nest.Connect(inputNeurons,spikedetectors)


### Training ###
trainingImages = sorted(glob.iglob("./images/training/*.png"))

img_fns = []
for img_fn in trainingImages:
    img_fns.append(img_fn)

if len(img_fns)!=nStim*nTrans:
    print 'Error: the number of images files does not match',len(img_fns);
    sys.exit(1)


for img_fn in img_fns:
    print img_fn
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
    
    

    plt.figure(1)
    plt.subplot(5,3,1);
    plt.imshow(img,interpolation='none');
    
    eccentricity = 30/10;
    
    rCounter = 1;
    for r in res_norm:   
        index = 0;  
        for neuron in inputNeurons:
            y_index = math.floor(index/layerDim);
            x_index = index%layerDim;
            #print mean(r[y_index][x_index]);
            nest.SetStatus([neuron], {"V_m": Vrest+(Vth-Vrest)*numpy.random.rand()})
            nest.SetStatus([neuron], {"I_e": 180.0+10*mean(r[y_index*3][x_index*3])}) #to-do: TOBE fixed
            index+=1;
            
        
        nest.SetStatus(spikedetectors, [{"n_events": 0}]);
        nest.Simulate(simulationTime)
        
        plt.subplot(5,3,rCounter*3+1);
        plt.imshow(res[rCounter-1],interpolation='none');
        
        dSD =nest.GetStatus(spikedetectors,keys='events')[0]
        evs = dSD["senders"]
        ts = dSD["times"]
        
        plt.subplot(5,3,rCounter*3+2);
        plt.plot(ts, evs, ".")
        
        res_FRMap = numpy.zeros((layerDim, layerDim));
        for index in evs:
            y_index = math.floor((index-1)/layerDim);
            x_index = (index-1)%layerDim;
            res_FRMap[y_index][x_index]+=1;
        
        res_FRMap*=(1000/simulationTime)
        
        plt.subplot(5,3,rCounter*3+3);
        plt.imshow(res_FRMap,interpolation='none');
        plt.colorbar();
        

        
        rCounter+=1;
    plt.show()
            
    

    
    #To-DO: trainNetworkWith
    
#Testing
