from Parameters import *
from imageImport import *
import sys
import numpy,pylab as plt , glob
import nest
import os

### Constructing Network ###
nest.Models()
ndict = {"I_e": I_e, "tau_m": tau_m}
#ndict = {"C_m":C_m, 'tau_m': tau_m, 't_ref': t_ref, 'E_L': E_L, 'V_th': V_th, 'V_reset': V_reset};
nest.SetDefaults("iaf_neuron", ndict)
inputNeurons = nest.Create("iaf_neuron", layerDim*layerDim*len(thetaList));


#nest.SetDefaults("iaf_psc_delta", ndict)
#inputNeurons = nest.Create('iaf_psc_alpha',layerDim*layerDim*len(thetaList))





spikedetectors = nest.Create("spike_detector", layerDim*layerDim*len(thetaList), params={"withgid": True, "withtime": True})
nest.Connect(inputNeurons,spikedetectors)


### Training ###
#trainingImages = sorted(glob.iglob("/images/training/*.png"))
trainingImages = sorted(glob.iglob(os.path.split(os.path.realpath(__file__))[0] + "/images/training/*.png"))

img_fns = []
for img_fn in trainingImages:
    img_fns.append(img_fn)

if len(img_fns)!=nStim*nTrans:
    print 'Error: the number of images files does not match',len(img_fns);
    sys.exit(1)

index_img=0;
for img_fn in img_fns:
    print img_fn
    #load gabor filtered ImageSurface
    img = cv2.imread(img_fn)
    if img is None:
        print 'Failed to load image file:', img_fn
        sys.exit(1)
    
    filters = build_filters()
    res = process(img, filters)
    
#     for r in res:
#         r= cv2.resize(r,(layerDim,layerDim),interpolation = cv2.INTER_NEAREST);
#     
    
    #convert from gabor filtered inputs to spikes
    #normalize
    res_norm=res/numpy.max(res);
    res_norm=1-res_norm;
    
    

    plt.figure(1)
    #plot input Image
    plt.subplot(5,3,1);
    plt.imshow(img,interpolation='none');
    plt.title('Input')
    
#     eccentricity = 30/10;
#     
#     plt.imshow(cv2.resize(res[2],(10,10),interpolation = cv2.INTER_NEAREST ),interpolation='none');
#     plt.show()
    
    
    #rCounter = 1;
    for index_filter in range(0,len(thetaList)):
        r = res_norm[index_filter]
        #index = 0;
        for index_cell in range(0,layerDim*layerDim-1):
        #for neuron in inputNeurons:
            y_index = math.floor(index_cell/layerDim);
            x_index = index_cell%layerDim;
            #print mean(r[y_index][x_index]);
            neuron = inputNeurons[index_filter*layerDim*layerDim+index_cell]
            nest.SetStatus([neuron], {"V_m": E_L+(V_th-E_L)*numpy.random.rand()})
            nest.SetStatus([neuron], {"I_e": I_e+80*mean(r[y_index][x_index])}) #TO-DO: TOBE fixed
            #print index_filter*len(thetaList)+index_cell
        
            #index+=1;
            
        
    nest.SetStatus(spikedetectors, [{"n_events": 0}]);
    nest.Simulate(simulationTime)

    dSD =nest.GetStatus(spikedetectors,keys='events')[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    
    for index_filter in range(0,len(res_norm)):    
        ax = plt.subplot(5,3,(index_filter+1)*3+1);
        plt.imshow(res[index_filter],interpolation='none');
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)
        plt.ylabel('Filter '+str(index_filter))
        
        evsIndex = [];
        res_FRMap = numpy.zeros((layerDim, layerDim));
        index_evs = 0;
        for index in evs:
            if layerDim*layerDim*index_filter <index and index<=layerDim*layerDim*(index_filter+1):
                evsIndex.append(index_evs);
                index_tmp = (index-1)%(layerDim*layerDim)+1
                y_index = math.floor((index_tmp-1)/layerDim);
                x_index = (index_tmp-1)%layerDim;
                res_FRMap[y_index][x_index]+=1;
            index_evs+=1;
        
        res_FRMap*=(1000/simulationTime)
        
        
        #plot FR map
        plt.subplot(5,3,(index_filter+1)*3+3);
        if(index_filter==0):
            plt.title('Firing Rate Map')
        plt.imshow(res_FRMap,interpolation='none');
        plt.colorbar();
        
        
        #plot spike raster
        ax=plt.subplot(5,3,(index_filter+1)*3+2);
        if(index_filter==0):
            plt.title('Raster Plot')
        plt.plot(ts[evsIndex], evs[evsIndex],'.')
        ax.get_yaxis().set_visible(False)
        ax.set_xlim([index_img*simulationTime, (index_img+1)*simulationTime])
        


    plt.show()
    index_img+=1;
            
    

    
    #To-DO: trainNetworkWith
    
#Testing
