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
nest.SetDefaults("stdp_synapse",stdp_dict);


plotGabor = 0;
plotLayer = 0; #0:each images, 1: at end

#Creating Gabor input layer
layerG = []
for theta in range(0,len(thetaList)):
    layerG.append(topp.CreateLayer(layerGDict))

#Creating ExcitLayers
layers = []
for layer in range(0,nLayers):
    layers.append(topp.CreateLayer(layersDict[layer]))

#Creating InhibLayers
inhibLayers = [] 
for layer in range(0,nLayers):
    inhibLayers.append(topp.CreateLayer(inhibLayersDict[layer]))

#Connecting neurons in Gabor input layer and neurons in the first ExcitLayer
for theta in range(0,len(thetaList)): 
    topp.ConnectLayers(layerG[theta], layers[0], connGDict)

#Connecting neurons between layers
for layer in range(0,nLayers-1): 
    topp.ConnectLayers(layers[layer], layers[layer+1], connForwardDict)#feed forward connections
    topp.ConnectLayers(layers[layer+1], layers[layer], connBackwardDict)#feedback connections

#Connecting neurons within layers    
for layer in range(0,nLayers): 
    topp.ConnectLayers(layers[layer], inhibLayers[layer], connExInhibDict)#E->I connections within layer 
    topp.ConnectLayers(inhibLayers[layer], layers[layer], connInhibExDict)#I->E connections within layer
    topp.ConnectLayers(inhibLayers[layer], inhibLayers[layer], connInhibRecDict)#I->I connections within layer

spkdetG = nest.Create("spike_detector", len(thetaList), params={"withgid": True, "withtime": True})
for theta in range(0,len(thetaList)):
    nest.Connect(nest.GetNodes(layerG[theta])[0],[spkdetG[theta]],"all_to_all");

spkdetLayers = []
for layer in range(0,nLayers):
    spkdetLayers.append(nest.Create("spike_detector", 1, params={"withgid": True, "withtime": True}))

for layer in range(0,nLayers):
    nest.Connect(nest.GetNodes(layers[layer])[0], [spkdetLayers[layer][0]], "all_to_all")

spkdetInhibLayers = []
for layer in range(0,nLayers):
    spkdetInhibLayers.append(nest.Create("spike_detector", 1, params={"withgid": True, "withtime": True}))

for layer in range(0,nLayers):
    nest.Connect(nest.GetNodes(inhibLayers[layer])[0], [spkdetInhibLayers[layer][0]], "all_to_all")

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
    
    
    for index_filter in range(0,len(thetaList)):
        nodesG=nest.GetNodes(layerG[index_filter])[0]

        r = res_norm[index_filter]
        for index_cell in range(0,layerGDim*layerGDim):
            y_index = int(math.floor(index_cell/layerGDim));
            x_index = index_cell%layerGDim;
            neuron = nodesG[index_cell]
            nest.SetStatus([neuron], {"V_m": E_L+(V_th-E_L)*numpy.random.rand()})
            nest.SetStatus([neuron], {"I_e": I_e+80*mean(r[y_index][x_index])}) #TO-DO: TOBE fixed

    if(plotGabor or plotLayer==0):
        nest.SetStatus(spkdetG, [{"n_events": 0}]);
        for layer in range(0,nLayers):
            nest.SetStatus(spkdetLayers[layer], [{"n_events": 0}])
            nest.SetStatus(spkdetInhibLayers[layer], [{"n_events": 0}])
    
    
    nest.Simulate(simulationTime)
    
    if(plotGabor):
        plt.figure(1)
        #plot input Image
        plt.subplot(5,3,1);
        plt.imshow(img,interpolation='none');
        plt.title('Input')
    
        for index_filter in range(0,len(res_norm)):    
            ax = plt.subplot(5,3,(index_filter+1)*3+1);
            plt.imshow(res[index_filter],interpolation='none');
            #ax.get_xaxis().set_visible(False)
            #ax.get_yaxis().set_visible(False)
            plt.ylabel('Filter '+str(index_filter))
            
            
            res_FRMap = numpy.zeros((layerGDim, layerGDim));
            
            dSD =nest.GetStatus(spkdetG,keys='events')[index_filter];
            evs = dSD["senders"]
            ts = dSD["times"]
                
            headNodeIndex = layerG[index_filter][0]
            for sender in evs:
                cell_index=sender-headNodeIndex-1
                y_index = math.floor((cell_index)/layerGDim);
                x_index = (cell_index)%layerGDim;
                res_FRMap[y_index][x_index]+=1*(1000/simulationTime);
         
            #plot spike raster
            ax=plt.subplot(5,3,(index_filter+1)*3+2);
            if(index_filter==0):
                plt.title('Raster Plot')
            plt.plot(ts, evs,'.')
            ax.get_yaxis().set_visible(False)
            ax.set_xlim([index_img*simulationTime, (index_img+1)*simulationTime])
    
            #plot FR map
            plt.subplot(5,3,(index_filter+1)*3+3);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(res_FRMap,interpolation='none',vmin=0, vmax=res_FRMap.max());
            plt.colorbar();
        
    if (plotLayer==0):
        plt.figure(2);    
        plt.subplot(nLayers+1,3,1);
        plt.imshow(img,interpolation='none');
        plt.title('Input')
        
        for layer in range(0,nLayers):
            
            ex_FRMap = numpy.zeros((layer1Dim, layer1Dim));
            dSD =nest.GetStatus(spkdetLayers[layer],keys='events')[0];
            evs = dSD["senders"]
            ts = dSD["times"]
                
            headNodeIndex = layers[layer][0]
            for sender in evs:
                cell_index=sender-headNodeIndex-1
                y_index = math.floor((cell_index)/layer1Dim);
                x_index = (cell_index)%layer1Dim;
                ex_FRMap[y_index][x_index]+=1*(1000/simulationTime);
         
            #plot spike raster
            ax=plt.subplot(nLayers+1,4,(nLayers-layer)*4+1);
            
            plt.title(layer)
            plt.plot(ts, evs,'.')
            #ax.get_yaxis().set_visible(False)
            ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
            ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);
            
            #plot FR map
            plt.subplot(nLayers+1,4,(nLayers-layer)*4+2);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(ex_FRMap,interpolation='none',vmin=0, vmax=ex_FRMap.max());
            plt.colorbar();
            
            
            
            
            inhib_FRMap = numpy.zeros((inhibLayer1Dim, inhibLayer1Dim));
            dSD =nest.GetStatus(spkdetInhibLayers[layer],keys='events')[0];
            evs = dSD["senders"]
            ts = dSD["times"]
                
            headNodeIndex = inhibLayers[layer][0]
            for sender in evs:
                cell_index=sender-headNodeIndex-1
                y_index = math.floor((cell_index)/inhibLayer1Dim);
                x_index = (cell_index)%inhibLayer1Dim;
                inhib_FRMap[y_index][x_index]+=1*(1000/simulationTime);
         
            #plot spike raster
            ax=plt.subplot(nLayers+1,4,(nLayers-layer)*4+3);
            
            plt.title(layer)
            plt.plot(ts, evs,'.')
            #ax.get_yaxis().set_visible(False)
            ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
            ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
            
            #plot FR map
            plt.subplot(nLayers+1,4,(nLayers-layer)*4+4);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(inhib_FRMap,interpolation='none',vmin=0, vmax=inhib_FRMap.max());
            plt.colorbar();
            

        plt.show()
    index_img+=1;
    
#To-DO: trainNetworkWith
if (plotLayer==1):
#         plt.figure(2);    
#         plt.subplot(nLayers+1,3,1);
#         plt.imshow(img,interpolation='none');
#         plt.title('Input')
    
    for layer in range(0,nLayers):
        headNodeIndex = layers[layer][0]
        dSD =nest.GetStatus(spkdetLayers[layer],keys='events')[0];
        evs = dSD["senders"]
        ts = dSD["times"]

        #plot spike raster
        ax=plt.subplot(nLayers,4,(nLayers-layer-1)*4+1);
        
        plt.title(layer)
        plt.plot(ts, evs,'.')
        #ax.get_yaxis().set_visible(False)
        #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
        ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);
        
        
        headNodeIndex = inhibLayers[layer][0]
        dSD =nest.GetStatus(spkdetInhibLayers[layer],keys='events')[0];
        evs = dSD["senders"]
        ts = dSD["times"]
            
        #plot spike raster
        ax=plt.subplot(nLayers,4,(nLayers-layer-1)*4+3);
        plt.title(layer)
        plt.plot(ts, evs,'.')
        #ax.get_yaxis().set_visible(False)
        #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
        ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
        
    plt.show()
    
#Testing
