from Parameters import *
from imageImport import *
import sys
import numpy,pylab as plt , glob
#import nest
import os
#from nest import raster_plot
#import nest.topology as topp

#http://brian2.readthedocs.org/en/latest/resources/tutorials/1-intro-to-brian-neurons.html

plotGabor = 1;
plotLayer = 0; #0:each images, 1: at end

net = Network(collect())

#Creating Gabor input layer
layerG = [];
for theta in range(0,len(thetaList)):
    #layerG.append(NeuronGroup(layerGDim*layerGDim, eqs, threshold='v>1', reset='v = 0'))
    layerG.append(PoissonGroup(layerGDim*layerGDim,numpy.random.rand(layerGDim*layerGDim)*Hz ))
net.add(layerG);

#Creating ExcitLayers
layers = []
for layer in range(0,nLayers):
    #layers.append(NeuronGroup(layerDim*layerDim, eqs, threshold='v>1', reset='v = 0'))
    layers.append(NeuronGroup(layerDim*layerDim, eqn_membran, threshold='v>Vt', reset='v = Vr', refractory=2*ms))
    layers[layer].v = 'Vr'# + rand() * (Vt - Vr)'
    layers[layer].ve = 0*mV
    layers[layer].vi = 0*mV

net.add(layers);

#Creating InhibLayers
inhibLayers = [] 
for layer in range(0,nLayers):
    #inhibLayers.append(NeuronGroup(inhibLayerDim*inhibLayerDim, eqs, threshold='v>1', reset='v = 0'))
    inhibLayers.append(NeuronGroup(layerDim*layerDim, eqn_membran, threshold='v>Vt', reset='v = Vr', refractory=2*ms))
    inhibLayers[layer].v = 'Vr'# + rand() * (Vt - Vr)'
    inhibLayers[layer].ve = 0*mV
    inhibLayers[layer].vi = 0*mV
net.add(inhibLayers);


#Connecting neurons in Gabor input layer and neurons in the first ExcitLayer
connGtoInput = []
for theta in range(0,len(thetaList)):
<<<<<<< local
    connGtoInput.append(Synapses(layerG[theta], layers[0], 'w:1', pre='ve += w*we'));
    connGtoInput[theta].w = 10;
    connGtoInput[theta].connect(True, p=0.2)
=======
    #connGtoInput.append(Synapses(layerG[theta], layers[0], 'w:1', pre='ve += w*we'));
    #connGtoInput[theta].w = 1;
    connGtoInput.append(Synapses(layerG[theta], layers[0], pre='ve += 10*we'));
    connGtoInput[theta].connect(True, p=0.1)
net.add(connGtoInput)
>>>>>>> other

#Connecting neurons between layers
connFeedForward = []
connBackProjection = []
for layer in range(0,nLayers-1): 
    connFeedForward.append(Synapses(layers[layer],layers[layer+1], 'w:1',pre='ve += 10*we'))
    #connFeedForward[layer].w = 1;
    connFeedForward[layer].connect(True, p=0.1);
    
    connBackProjection.append(Synapses(layers[layer+1],layers[layer], 'w:1',pre='ve += 10*we'))
    #connBackProjection[layer].w = 1;
    connBackProjection[layer].connect(True, p=0.1);
net.add(connFeedForward)
net.add(connBackProjection)
    
#Connecting neurons within layers    
connExIn = []
connInEx = []
connRecIn = []
connRecEx = []
for layer in range(0,nLayers): 
    connExIn.append(Synapses(layers[layer],inhibLayers[layer], 'w:1',pre='ve += 10*we'))
    #connExIn[layer].w = 1;
    connExIn[layer].connect(True, p=0.1)
    
    connInEx.append(Synapses(inhibLayers[layer], layers[layer],'w:1', pre='vi += 10*wi'))
    #connInEx[layer].w = 1;
    connInEx[layer].connect(True, p=0.1)
net.add(connExIn)
net.add(connInEx)



spikesG=[]
for theta in range(0,len(thetaList)):
    spikesG.append(SpikeMonitor(layerG[theta]))
net.add(spikesG)

# tmp = layerG[0];
# testSpikes = SpikeMonitor(tmp);

spkdetLayers = []
for layer in range(0,nLayers):
    spkdetLayers.append(SpikeMonitor(layers[layer]))
    #spkdetLayers.append(StateMonitor(layers[layer], 'v', record=True))
net.add(spkdetLayers)

spkdetInhibLayers = []
for layer in range(0,nLayers):
    spkdetInhibLayers.append(SpikeMonitor(inhibLayers[layer]))
net.add(spkdetInhibLayers)

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
                     
            #plot spike raster
            ax=plt.subplot(5,3,(index_filter+1)*3+2);
            if(index_filter==0):
                plt.title('Raster Plot')
            tmp = spikesG[index_filter];
            #plot(spikesG[index_filter].t/ms, spikesG[index_filter].i, '.')
            #plot(testSpikes.t/ms, testSpikes.i, '.')
            plot(tmp.t/ms, tmp.i, '.')
            #plot FR map
            plt.subplot(5,3,(index_filter+1)*3+3);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(res_FRMap,interpolation='none',vmin=0, vmax=res_FRMap.max());
            plt.colorbar();
        #plt.show();
        
        
    if (plotLayer==0):
        plt.figure(2);    
        plt.subplot(nLayers+1,3,1);
        plt.imshow(img,interpolation='none');
        plt.title('Input')
        
        for layer in range(0,nLayers):
            
            ex_FRMap = numpy.zeros((layerDim, layerDim));
            #plot spike raster
            ax=plt.subplot(nLayers+1,4,(nLayers-layer)*4+1);
            
            plt.title(layer)
            #plot(spkdetLayers[layer].t/ms, spkdetLayers[layer].v[0])
            plot(spkdetLayers[layer].t/ms, spkdetLayers[layer].i, '.')
            #ax.get_yaxis().set_visible(False)
            #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
            #ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);
            
            #plot FR map
            plt.subplot(nLayers+1,4,(nLayers-layer)*4+2);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(ex_FRMap,interpolation='none',vmin=0, vmax=ex_FRMap.max());
            plt.colorbar();
            
          
            #plot spike raster
            ax=plt.subplot(nLayers+1,4,(nLayers-layer)*4+3);
             
            plt.title(layer)
            plot(spkdetInhibLayers[layer].t/ms, spkdetInhibLayers[layer].i, '.')
#             plt.plot(ts, evs,'.')
#             #ax.get_yaxis().set_visible(False)
#             ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
#             ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
#             
#             #plot FR map
#             plt.subplot(nLayers+1,4,(nLayers-layer)*4+4);
#             if(index_filter==0):
#                 plt.title('Firing Rate Map')
#             plt.imshow(inhib_FRMap,interpolation='none',vmin=0, vmax=inhib_FRMap.max());
#             plt.colorbar();
            

        plt.show()
    index_img+=1;
    
#To-DO: trainNetworkWith
if (plotLayer==1):
#         plt.figure(2);    
#         plt.subplot(nLayers+1,3,1);
#         plt.imshow(img,interpolation='none');
#         plt.title('Input')
    
    for layer in range(0,nLayers):
        #headNodeIndex = layers[layer][0]

        #plot spike raster
        ax=plt.subplot(nLayers,4,(nLayers-layer-1)*4+1);
        
        plt.title(layer)
        raster_plot(spkdetLayers[layer])
        #ax.get_yaxis().set_visible(False)
        #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
        #ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);
        
        
        #dSD =nest.GetStatus(spkdetInhibLayers[layer],keys='events')[0];
            
        #plot spike raster
        ax=plt.subplot(nLayers,4,(nLayers-layer-1)*4+3);
        plt.title(layer)
        raster_plot(spkdetInhibLayers[layer])
        #ax.get_yaxis().set_visible(False)
        #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
        #ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
        
    plt.show()
    
#Testing
