import pylab as plt
import numpy
from Parameters import *
from CreateNetwork import *

plotGabor = 1;
plotLayer = 1; #1:each stimulus, 2: at end

def PlotGabor (img,res_norm,res, index_img):
    
    if(plotGabor):
        plt.figure(1 , figsize=(20,10))
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
            plt.ylim([0,layerGDim*layerGDim-1])
            plt.xlim([simulationTime*index_img,simulationTime*(index_img+1)])
        
        
            tmp2 = spikesG[index_filter].spike_trains();
            for row_tmp in range(layerGDim):
                for col_tmp in range(layerGDim):
                    index_tmp = row_tmp*layerGDim + col_tmp;
                    #print tmp2[index_tmp]
                    
                    if (len(tmp2[index_tmp])==0):
                        condition = tmp2[index_tmp]>simulationTime*index_img
                    else:
                        condition = tmp2[index_tmp]>simulationTime*index_img*ms
                
                    res_FRMap[row_tmp][col_tmp] = len(numpy.extract(condition,tmp2[index_tmp]));
        
            #plot FR map
            plt.subplot(5,3,(index_filter+1)*3+3);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(res_FRMap,interpolation='none',vmin=0, vmax=Rmax);
            plt.colorbar();
    #plt.show();
    return 0

def PlotLayer (img,res_norm,res, index_img, index_filter):
    if (plotLayer==1):
        plt.figure(2, figsize=(20,10) );
        plt.subplot(nLayers+1,3,1);
        plt.imshow(img,interpolation='none');
        plt.title('Input')
        
        for layer in range(0,nLayers):
            
            ex_FRMap = numpy.zeros((layerDim, layerDim));
            #plot spike raster
            ax=plt.subplot(nLayers+1,4,(nLayers-layer)*4+1);
            
            plt.title('Excitatory layer '+str(layer))
            #plot(spkdetLayers[layer].t/ms, spkdetLayers[layer].v[0])
            plot(spkdetLayers[layer].t/ms, spkdetLayers[layer].i, '.')
            plt.ylim([0,layerGDim*layerGDim-1])
            plt.xlim([simulationTime*index_img,simulationTime*(index_img+1)])
            #ax.get_yaxis().set_visible(False)
            #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
            #ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);
            
            #plot FR map
            ex_FRMap = numpy.zeros((layerDim, layerDim));
            tmp2 = spkdetInhibLayers[layer].spike_trains();
            for row_tmp in range(layerDim):
                for col_tmp in range(layerDim):
                    index_tmp = row_tmp*layerGDim + col_tmp;
                    
                    if (len(tmp2[index_tmp])==0):
                        condition = tmp2[index_tmp]>simulationTime*index_img
                    else:
                        condition = tmp2[index_tmp]>simulationTime*index_img*ms
                
                    ex_FRMap[row_tmp][col_tmp] = len(numpy.extract(condition,tmp2[index_tmp]));
            
            plt.subplot(nLayers+1,4,(nLayers-layer)*4+2);
            if(index_filter==0):
                plt.title('Firing Rate Map')
            plt.imshow(ex_FRMap,interpolation='none',vmin=0, vmax=ex_FRMap.max());
            plt.colorbar();
        
        
            #plot spike raster
            ax=plt.subplot(nLayers+1,4,(nLayers-layer)*4+3);
        
            plt.title('Inhibitory layer '+str(layer))
            plt.ylim([0,layerGDim*layerGDim-1])
            plt.xlim([simulationTime*index_img,simulationTime*(index_img+1)])
            plot(spkdetInhibLayers[layer].t/ms, spkdetInhibLayers[layer].i, '.')
            #             plt.plot(ts, evs,'.')
            #             #ax.get_yaxis().set_visible(False)
            #             ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
            #             ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
            #
            #plot FR map
        
            inhib_FRMap = numpy.zeros((layerDim, layerDim));
            tmp2 = spkdetInhibLayers[layer].spike_trains();
            for row_tmp in range(layerDim):
                for col_tmp in range(layerDim):
                    index_tmp = row_tmp*layerGDim + col_tmp;
                    inhib_FRMap[row_tmp][col_tmp] = len(tmp2[index_tmp]);
        
            
            plt.subplot(nLayers+1,4,(nLayers-layer)*4+4);
            
            plt.title('Firing Rate Map')
            plt.imshow(inhib_FRMap,interpolation='none',vmin=0, vmax=inhib_FRMap.max());
            plt.colorbar();
        
        plt.show()


    elif (plotLayer==2):
        
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
            pylab.ylim([0,layerGDim*layerGDim-1])
            pylab.xlim([simulationTime*index_img,simulationTime*index_img+1])
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

        
    return 0
