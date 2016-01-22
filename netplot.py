# date: 24/11/15
# author: Akihiro Eguchi
# description: a class to plot results

import pylab as plt

#from Parameters import *


class plotter(object):

    """a plotter for spiking visnet simulations"""

    def __init__(self, vnet,borrowed_globals):
        super(plotter, self).__init__()
        globals().update(borrowed_globals);

        self.vnet = vnet
        self.timeBegin = 0;

    def plotGaborInput(self, img, index_img, res, res_norm,simulationTime):

        self.figG = plt.figure(1 , figsize=(20, 10))
        # plot input Image
        plt.subplot(5, 3, 1)
        plt.imshow(img, cmap='gray', vmin=0, vmax=255, interpolation='none')
        plt.title('Input')


        
        ps = len(psiList);
        ss = len(lamdaList);
        ors = len(thetaList);

        for p in range(ps):
            for s in range(ss):
                for o in range(ors):
                    index_filter = p * (ss * ors) + s * ors + o;
                    ax = plt.subplot(5, 3, (index_filter + 1) * 3 + 1)
                    plt.imshow(res[p][s][o], cmap='jet', interpolation='none')
                    # ax.get_xaxis().set_visible(False)
                    # ax.get_yaxis().set_visible(False)
                    plt.ylabel('Filter ' + str(index_filter))
        
                    res_FRMap = np.zeros((layerGDim, layerGDim))
        
                    # plot spike raster
                    ax = plt.subplot(5, 3, (index_filter + 1) * 3 + 2)
                    if(index_filter == 0):
                        plt.title('Raster Plot')
                    tmp = self.vnet.spikesG[index_filter]
                    # plot(self.vnet.spikesG[index_filter].t/ms, self.vnet.spikesG[index_filter].i, '.')
                    # plot(testSpikes.t/ms, testSpikes.i, '.')
                    plt.plot(tmp.t / ms, tmp.i, '.')
                    plt.ylim([0, layerGDim * layerGDim - 1])
                    plt.xlim([self.timeBegin, self.timeBegin+simulationTime])
        
                    tmp2 = self.vnet.spikesG[index_filter].spike_trains()
                    for row_tmp in range(layerGDim):
                        for col_tmp in range(layerGDim):
                            index_tmp = row_tmp * layerGDim + col_tmp
                            if (len(tmp2[index_tmp]) == 0):
                                condition = tmp2[index_tmp] > self.timeBegin
                            else:
                                condition = tmp2[index_tmp] > self.timeBegin * ms
                            res_FRMap[row_tmp][col_tmp] = len(np.extract(condition, tmp2[index_tmp]));
        
                    res_FRMap = res_FRMap/(float(simulationTime)/1000);
                    
                    # plot FR map
                    plt.subplot(5, 3, (index_filter + 1) * 3 + 3)
                    if(index_filter == 0):
                        plt.title('Firing Rate Map')
                    plt.imshow(res_FRMap, cmap='jet', interpolation='none', vmin=0, vmax=Rmax)
                    plt.colorbar()

    def plotLayers(self, img, index_img, FRrecTmp, simulationTime):
        self.figL = plt.figure(2, figsize=(20, 10));
        plt.title('Input')
        plt.subplot(nLayers + 1, 3, 1)
        plt.imshow(img, cmap='gray', vmin=0, vmax=255, interpolation='none')
        
        # plot binding layer
        # plot spike raster of binding layer
        plt.subplot(nLayers + 1, 3, 2)
        plt.title('Binding layer')
        plt.plot(
            self.vnet.spkdetBindingLayer.t / ms,
            self.vnet.spkdetBindingLayer.i, '.')
        plt.ylim([0, layerDim * layerDim - 1])
        plt.xlim([self.timeBegin, self.timeBegin+simulationTime])
        
        # plot FR map of bindingLayer
        bind_FRMap = np.zeros((layerDim, layerDim))
        tmp2 = self.vnet.spkdetBindingLayer.spike_trains()
        for row_tmp in range(layerDim):
            for col_tmp in range(layerDim):
                index_tmp = row_tmp * layerDim + col_tmp
                if (len(tmp2[index_tmp]) == 0):
                    condition = tmp2[index_tmp] > self.timeBegin
                else:
                    condition = tmp2[index_tmp] > self.timeBegin * ms
                bind_FRMap[row_tmp, col_tmp] = len(np.extract(condition, tmp2[index_tmp]));
        bind_FRMap = bind_FRMap/(float(simulationTime)/1000);
        
        plt.subplot(nLayers + 1, 3, 3)
        plt.imshow(
            bind_FRMap, cmap='jet', interpolation='none', vmin=0, vmax=bind_FRMap.max())
        plt.colorbar()
        FRrecTmp[0] = bind_FRMap;
        
        

        for layer in range(0, nLayers):
            # plot spike raster of excitatory
            ax = plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 1)
            plt.title('Excitatory layer ' + str(layer))
            plt.plot(
                self.vnet.spkdetLayers[layer].t / ms,
                self.vnet.spkdetLayers[layer].i, '.')
            plt.ylim([0, layerDim * layerDim - 1])
            plt.xlim([self.timeBegin, self.timeBegin+simulationTime])
            
            # plot FR map of excitatory
            ex_FRMap = np.zeros((layerDim, layerDim))
            tmp2 = self.vnet.spkdetLayers[layer].spike_trains()
            for row_tmp in range(layerDim):
                for col_tmp in range(layerDim):
                    index_tmp = row_tmp * layerDim + col_tmp
                    if (len(tmp2[index_tmp]) == 0):
                        condition = tmp2[index_tmp] > self.timeBegin
                    else:
                        condition = tmp2[index_tmp] > self.timeBegin * ms
                    ex_FRMap[row_tmp, col_tmp] = len(np.extract(condition, tmp2[index_tmp]));
            ex_FRMap = ex_FRMap/(float(simulationTime)/1000);
            
            plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 2)
            plt.imshow(
                ex_FRMap, cmap='jet', interpolation='none', vmin=0, vmax=ex_FRMap.max())
            plt.colorbar()
            FRrecTmp[layer + 1] = ex_FRMap;
            
            
            # plot spike raster of inhibitory
            ax = plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 3)

            plt.title('Inhibitory layer ' + str(layer))
            plt.ylim([0, inhibLayerDim * inhibLayerDim - 1])
            plt.xlim([self.timeBegin, self.timeBegin+simulationTime])
            plt.plot(
                self.vnet.spkdetInhibLayers[layer].t / ms,
                self.vnet.spkdetInhibLayers[layer].i, '.')

            # plot FR map of inhibitory
            inhib_FRMap = np.zeros((inhibLayerDim, inhibLayerDim))
            tmp2 = self.vnet.spkdetInhibLayers[layer].spike_trains()
            for row_tmp in range(inhibLayerDim):
                for col_tmp in range(inhibLayerDim):
                    index_tmp = row_tmp * inhibLayerDim + col_tmp
                    inhib_FRMap[row_tmp, col_tmp] = len(tmp2[index_tmp])
            inhib_FRMap = inhib_FRMap/(float(simulationTime)/1000);

            plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 4)
            plt.imshow(
                inhib_FRMap, cmap='jet', interpolation='none', vmin=0,
                vmax=inhib_FRMap.max())
            plt.title('Firing Rate Map')
            plt.colorbar()

        
        self.timeBegin += simulationTime;
        # plt.show()
    
    def saveFigs(self,plotActivities,plotGabor,svname):
        if plotActivities:
            self.figL.savefig(svname+"_l.png");
        if plotGabor:
            self.figG.savefig(svname+"_g.png");
        
        
        
