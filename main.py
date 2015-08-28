from Parameters import *
import imageImport as imimp
import cv2
import sys
import numpy,pylab as plt , glob
#import nest
import os
import netarch
#from nest import raster_plot
#import nest.topology as topp
plotGabor = 1;
plotLayer = 0; #0:each images, 1: at end

#http://brian2.readthedocs.org/en/latest/resources/tutorials/1-intro-to-brian-neurons.html

# build visnet
vnet = netarch.visnet()

### Training ###
#trainingImages = sorted(glob.iglob("/images/training/*.png"))
img_fns = sorted(glob.iglob(os.path.split(os.path.realpath(__file__))[0] + "/images/training/*.png"))

# img_fns = []
# for img_fn in trainingImages:
#     img_fns.append(img_fn)

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
        r = numpy.reshape(numpy.mean(res_norm[index_filter],axis=2),(layerGDim*layerGDim));
        #print r
        vnet.layerG[index_filter].rates= r * Rmax;    #To be fixed
        #print vnet.layerG[index_filter].rates
        
    
    
    
    vnet.net.run(simulationTime*ms)
    
    
    
    
    
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
            tmp = vnet.spikesG[index_filter];
            #plot(vnet.spikesG[index_filter].t/ms, vnet.spikesG[index_filter].i, '.')
            #plot(testSpikes.t/ms, testSpikes.i, '.')
            plt.plot(tmp.t/ms, tmp.i, '.')
            
            
            tmp2 = vnet.spikesG[index_filter].spike_trains();
            for row_tmp in range(layerGDim):
                for col_tmp in range(layerGDim):
                    index_tmp = row_tmp*layerGDim + col_tmp;
                    res_FRMap[row_tmp][col_tmp] = len(tmp2[index_tmp]);
            
            #plot FR map
            plt.subplot(5,3,(index_filter+1)*3+3);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(res_FRMap,interpolation='none',vmin=0, vmax=Rmax);
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
            #plot(vnet.spkdetLayers[layer].t/ms, vnet.spkdetLayers[layer].v[0])
            plt.plot(vnet.spkdetLayers[layer].t/ms, vnet.spkdetLayers[layer].i, '.')
            #ax.get_yaxis().set_visible(False)
            #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
            #ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);
            
            #plot FR map
            ex_FRMap = numpy.zeros((layerDim, layerDim));
            tmp2 = vnet.spkdetInhibLayers[layer].spike_trains();
            for row_tmp in range(layerDim):
                for col_tmp in range(layerDim):
                    index_tmp = row_tmp*layerGDim + col_tmp;
                    ex_FRMap[row_tmp][col_tmp] = len(tmp2[index_tmp]);
                    
            plt.subplot(nLayers+1,4,(nLayers-layer)*4+2);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(ex_FRMap,interpolation='none',vmin=0, vmax=ex_FRMap.max());
            plt.colorbar();
            
          
            #plot spike raster
            ax=plt.subplot(nLayers+1,4,(nLayers-layer)*4+3);
             
            plt.title(layer)
            plt.plot(vnet.spkdetInhibLayers[layer].t/ms, vnet.spkdetInhibLayers[layer].i, '.')
#             plt.plot(ts, evs,'.')
#             #ax.get_yaxis().set_visible(False)
#             ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
#             ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
#             
            #plot FR map
            
            inhib_FRMap = numpy.zeros((layerDim, layerDim));
            tmp2 = vnet.spkdetInhibLayers[layer].spike_trains();
            for row_tmp in range(layerDim):
                for col_tmp in range(layerDim):
                    index_tmp = row_tmp*layerGDim + col_tmp;
                    inhib_FRMap[row_tmp][col_tmp] = len(tmp2[index_tmp]);
            
            
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
        #headNodeIndex = vnet.layers[layer][0]

        #plot spike raster
        ax=plt.subplot(nLayers,4,(nLayers-layer-1)*4+1);
        
        plt.title(layer)
        raster_plot(vnet.spkdetLayers[layer])
        #ax.get_yaxis().set_visible(False)
        #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
        #ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);
        
        
        #dSD =nest.GetStatus(vnet.spkdetInhibLayers[layer],keys='events')[0];
            
        #plot spike raster
        ax=plt.subplot(nLayers,4,(nLayers-layer-1)*4+3);
        plt.title(layer)
        raster_plot(vnet.spkdetInhibLayers[layer])
        #ax.get_yaxis().set_visible(False)
        #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
        #ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
        
    plt.show()
    
#Testing
