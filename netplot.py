
import pylab as plt

from Parameters import *


class plotter(object):

    """a plotter for spiking visnet simulations"""

    def __init__(self, vnet):
        super(plotter, self).__init__()

        self.vnet = vnet

    def plotGaborInput(self, img, index_img, res, res_norm):

        plt.figure(1 , figsize=(20, 10))
        # plot input Image
        plt.subplot(5, 3, 1)
        plt.imshow(img, interpolation='none')
        plt.title('Input')

        for index_filter in range(0, len(res_norm)):
            ax = plt.subplot(5, 3, (index_filter + 1) * 3 + 1)
            plt.imshow(res[index_filter], interpolation='none')
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
            plt.xlim([simulationTime * index_img, simulationTime * (index_img + 1)])

            tmp2 = self.vnet.spikesG[index_filter].spike_trains()
            for row_tmp in range(layerGDim):
                for col_tmp in range(layerGDim):
                    index_tmp = row_tmp * layerGDim + col_tmp
                    if (len(tmp2[index_tmp]) == 0):
                        condition = tmp2[index_tmp] > simulationTime * index_img
                    else:
                        condition = tmp2[index_tmp] > simulationTime * index_img * ms
                    res_FRMap[row_tmp][col_tmp] = len(np.extract(condition, tmp2[index_tmp]));

            # plot FR map
            plt.subplot(5, 3, (index_filter + 1) * 3 + 3)
            if(index_filter == 0):
                plt.title('Firing Rate Map')
            plt.imshow(res_FRMap, interpolation='none', vmin=0, vmax=Rmax)
            plt.colorbar()

    def plotLayers(self, img, index_img, FRrecTmp):

        # if (plotLayer == 0):
            # plot input image
            plt.figure(2, figsize=(20, 10));
            plt.title('Input')
            plt.subplot(nLayers + 1, 3, 1)
            plt.imshow(img, interpolation='none')
            
            # plot binding layer
            # plot spike raster of binding layer
            plt.subplot(nLayers + 1, 3, 2)
            plt.title('Binding layer')
            plt.plot(
                self.vnet.spkdetBindingLayer.t / ms,
                self.vnet.spkdetBindingLayer.i, '.')
            plt.ylim([0, layerGDim * layerGDim - 1])
            plt.xlim([simulationTime * index_img, simulationTime * (index_img + 1)])
            
            # plot FR map of bindingLayer
            bind_FRMap = np.zeros((layerDim, layerDim))
            tmp2 = self.vnet.spkdetBindingLayer.spike_trains()
            for row_tmp in range(layerDim):
                for col_tmp in range(layerDim):
                    index_tmp = row_tmp * layerGDim + col_tmp
                    if (len(tmp2[index_tmp]) == 0):
                        condition = tmp2[index_tmp] > simulationTime * index_img
                    else:
                        condition = tmp2[index_tmp] > simulationTime * index_img * ms
                    bind_FRMap[row_tmp,col_tmp] = len(np.extract(condition, tmp2[index_tmp]));
            plt.subplot(nLayers + 1, 3, 3)
            plt.imshow(
                bind_FRMap, interpolation='none', vmin=0, vmax=bind_FRMap.max())
            plt.colorbar()
            FRrecTmp[0]=bind_FRMap;
            
            

            for layer in range(0, nLayers):
                # plot spike raster of excitatory
                ax = plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 1)
                plt.title('Excitatory layer ' + str(layer))
                plt.plot(
                    self.vnet.spkdetLayers[layer].t / ms,
                    self.vnet.spkdetLayers[layer].i, '.')
                plt.ylim([0, layerGDim * layerGDim - 1])
                plt.xlim([simulationTime * index_img, simulationTime * (index_img + 1)])
                
                # plot FR map of excitatory
                ex_FRMap = np.zeros((layerDim, layerDim))
                tmp2 = self.vnet.spkdetLayers[layer].spike_trains()
                for row_tmp in range(layerDim):
                    for col_tmp in range(layerDim):
                        index_tmp = row_tmp * layerGDim + col_tmp
                        if (len(tmp2[index_tmp]) == 0):
                            condition = tmp2[index_tmp] > simulationTime * index_img
                        else:
                            condition = tmp2[index_tmp] > simulationTime * index_img * ms
                        ex_FRMap[row_tmp,col_tmp] = len(np.extract(condition, tmp2[index_tmp]));
                
                plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 2)
                plt.imshow(
                    ex_FRMap, interpolation='none', vmin=0, vmax=ex_FRMap.max())
                plt.colorbar()
                FRrecTmp[layer+1]=ex_FRMap;
                
                
                # plot spike raster of inhibitory
                ax = plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 3)

                plt.title('Inhibitory layer ' + str(layer))
                plt.ylim([0, layerGDim * layerGDim - 1])
                plt.xlim([simulationTime * index_img, simulationTime * (index_img + 1)])
                plt.plot(
                    self.vnet.spkdetInhibLayers[layer].t / ms,
                    self.vnet.spkdetInhibLayers[layer].i, '.')

                # plot FR map of inhibitory
                inhib_FRMap = np.zeros((layerDim, layerDim))
                tmp2 = self.vnet.spkdetInhibLayers[layer].spike_trains()
                for row_tmp in range(layerDim):
                    for col_tmp in range(layerDim):
                        index_tmp = row_tmp * layerGDim + col_tmp
                        inhib_FRMap[row_tmp,col_tmp] = len(tmp2[index_tmp])

                plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 4)
                plt.imshow(
                    inhib_FRMap, interpolation='none', vmin=0,
                    vmax=inhib_FRMap.max())
                plt.title('Firing Rate Map')
                plt.colorbar()

            #plt.show()
