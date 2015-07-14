from Parameters import *
from imageImport import *
import sys
import numpy,pylab as plt , glob
import nest
import os
import nest.topology as topp

### Constructing Network ###
nest.Models()
#ndict = {"C_m":C_m, 'tau_m': tau_m, 't_ref': t_ref, 'E_L': E_L, 'V_th': V_th, 'V_reset': V_reset};
nest.SetDefaults("iaf_neuron", ndict)

#layerG = nest.Create("iaf_neuron", layerDim*layerDim*len(thetaList));


layerG = topp.CreateLayer(layerGDict)
layer1 = topp.CreateLayer(layer1Dict)
layer2 = topp.CreateLayer(layer2Dict)

topp.ConnectLayers(layerG, layer1, connGDict)
topp.ConnectLayers(layer1, layer2, conn1Dict)
topp.ConnectLayers(layer2, layer1, conn2Dict)


# nest.SetDefaults("spike_detector", spkdet1Dict)
# spkdet1 = topp.CreateLayer(layer1SpkDict)
# nest.SetDefaults("spike_detector", spkdet2Dict)
# spkdet2 = topp.CreateLayer(layer2SpkDict)
# nest.SetDefaults("spike_detector", spkdetGDict)
# spkdetG = topp.CreateLayer(layerGSpkDict)

# topp.ConnectLayers(layerG, spkdetG, connSpkDict)
# topp.ConnectLayers(layer1, spkdet1, connSpkDict)
# topp.ConnectLayers(layer2, spkdet2, connSpkDict)

spikedetectors = nest.Create("spike_detector", layerGDim*layerGDim*len(thetaList), params={"withgid": True, "withtime": True})
nest.Connect(nest.GetNodes(layerG)[0],spikedetectors);


#nest.SetDefaults("iaf_psc_delta", ndict)
#inputNeurons = nest.Create('iaf_psc_alpha',layerDim*layerDim*len(thetaList))





#spikedetectors = nest.Create("spike_detector", layerDim*layerDim*len(thetaList), params={"withgid": True, "withtime": True})
#nest.Connect(inputNeurons,spikedetectors)


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
    #print img_fn
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
    
    nodesG=nest.GetNodes(layerG)[0]
    


    #rCounter = 1;
    for index_filter in range(0,len(thetaList)-1):
        r = res_norm[index_filter]
        #index = 0;
        for index_cell in range(0,layerGDim*layerGDim-1):
        #for neuron in inputNeurons:
            y_index = int(math.floor(index_cell/layerGDim));
            x_index = index_cell%layerGDim;
            #print mean(r[y_index][x_index]);
            neuron = nodesG[index_filter*layerGDim*layerGDim+index_cell]
            nest.SetStatus([neuron], {"V_m": E_L+(V_th-E_L)*numpy.random.rand()})
            nest.SetStatus([neuron], {"I_e": I_e+80*mean(r[y_index][x_index])}) #TO-DO: TOBE fixed
            #print index_filter*len(thetaList)+index_cell
        
            #index+=1;
            
    #nodesSpkG=nest.GetNodes(spkdetG)[0]
    #nest.SetStatus(nodesSpkG, {"n_events": 0});
    nest.SetStatus(spikedetectors, [{"n_events": 0}]);
    
    nest.Simulate(simulationTime)

    dSD =nest.GetStatus(spikedetectors,keys='events')[0]
    #dSD =nest.GetStatus(nodesSpkG,keys='events')[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    
    
    
    
    for index_filter in range(0,len(res_norm)):    
        ax = plt.subplot(5,3,(index_filter+1)*3+1);
        plt.imshow(res[index_filter],interpolation='none');
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)
        plt.ylabel('Filter '+str(index_filter))
        
        evsIndex = [];
        res_FRMap = numpy.zeros((layerGDim, layerGDim));
        index_evs = 0;
        for index in evs:
            if layerGDim*layerGDim*index_filter <index and index<=layerGDim*layerGDim*(index_filter+1):
                evsIndex.append(index_evs);
                index_tmp = (index-1)%(layerGDim*layerGDim)+1
                y_index = math.floor((index_tmp-1)/layerGDim);
                x_index = (index_tmp-1)%layerGDim;
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
