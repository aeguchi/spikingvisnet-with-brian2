
import pylab as plt
from Parameters import *


class plotter(object):

    """a plotter for spiking visnet simulations"""

    def __init__(self, vnet):
        super(plotter, self).__init__()

        self.vnet = vnet

    def plotGaborInput(self, img, res, res_norm):

        plt.figure(1)
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
            #plot(self.vnet.spikesG[index_filter].t/ms, self.vnet.spikesG[index_filter].i, '.')
            #plot(testSpikes.t/ms, testSpikes.i, '.')
            plt.plot(tmp.t / ms, tmp.i, '.')

            tmp2 = self.vnet.spikesG[index_filter].spike_trains()
            for row_tmp in range(layerGDim):
                for col_tmp in range(layerGDim):
                    index_tmp = row_tmp * layerGDim + col_tmp
                    res_FRMap[row_tmp][col_tmp] = len(tmp2[index_tmp])

            # plot FR map
            plt.subplot(5, 3, (index_filter + 1) * 3 + 3)
            if(index_filter == 0):
                plt.title('Firing Rate Map')
            plt.imshow(res_FRMap, interpolation='none', vmin=0, vmax=Rmax)
            plt.colorbar()

    def plotLayers(self, img, plotLayer=0):

        if (plotLayer == 0):

            plt.figure(2)
            plt.subplot(nLayers + 1, 3, 1)
            plt.imshow(img, interpolation='none')
            plt.title('Input')

            for layer in range(0, nLayers):

                ex_FRMap = np.zeros((layerDim, layerDim))
                # plot spike raster
                ax = plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 1)

                plt.title(layer)
                # plot(self.vnet.spkdetLayers[layer].t/ms,
                #      self.vnet.spkdetLayers[layer].v[0])
                plt.plot(
                    self.vnet.spkdetLayers[layer].t / ms,
                    self.vnet.spkdetLayers[layer].i, '.')
                # ax.get_yaxis().set_visible(False)
                # ax.set_xlim([(index_img)*simulationTime,
                #              (index_img+1)*simulationTime])
                # ax.set_ylim([headNodeIndex+1,
                #             headNodeIndex+(layer1Dim*layer1Dim)]);

                # plot FR map
                ex_FRMap = np.zeros((layerDim, layerDim))
                tmp2 = self.vnet.spkdetInhibLayers[layer].spike_trains()
                for row_tmp in range(layerDim):
                    for col_tmp in range(layerDim):
                        index_tmp = row_tmp * layerGDim + col_tmp
                        ex_FRMap[row_tmp][col_tmp] = len(tmp2[index_tmp])

                plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 2)
                # if(index_filter == 0):
                #     plt.title('Firing Rate Map')
                plt.imshow(
                    ex_FRMap, interpolation='none', vmin=0, vmax=ex_FRMap.max())
                plt.colorbar()

                # plot spike raster
                ax = plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 3)

                plt.title(layer)
                plt.plot(
                    self.vnet.spkdetInhibLayers[layer].t / ms,
                    self.vnet.spkdetInhibLayers[layer].i, '.')

    #             plt.plot(ts, evs,'.')
    # ax.get_yaxis().set_visible(False)
    #             ax.set_xlim([(index_img)*simulationTime,
    #                          (index_img+1)*simulationTime])
    #             ax.set_ylim([headNodeIndex+1,
    #                          headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
    #
                # plot FR map

                inhib_FRMap = np.zeros((layerDim, layerDim))
                tmp2 = self.vnet.spkdetInhibLayers[layer].spike_trains()
                for row_tmp in range(layerDim):
                    for col_tmp in range(layerDim):
                        index_tmp = row_tmp * layerGDim + col_tmp
                        inhib_FRMap[row_tmp][col_tmp] = len(tmp2[index_tmp])

                plt.subplot(nLayers + 1, 4, (nLayers - layer) * 4 + 4)
                # if(index_filter == 0):
                #     plt.title('Firing Rate Map')
                plt.imshow(
                    inhib_FRMap, interpolation='none', vmin=0,
                    vmax=inhib_FRMap.max())
                plt.colorbar()

            plt.show()

        # To-DO: trainNetworkWith
        if (plotLayer == 1):
        #         plt.figure(2);
        #         plt.subplot(nLayers+1,3,1);
        #         plt.imshow(img,interpolation='none');
        #         plt.title('Input')

            for layer in range(0, nLayers):
                #headNodeIndex = self.vnet.layers[layer][0]

                # plot spike raster
                ax = plt.subplot(nLayers, 4, (nLayers - layer - 1) * 4 + 1)

                plt.title(layer)
                raster_plot(self.vnet.spkdetLayers[layer])
                # ax.get_yaxis().set_visible(False)
                #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
                #ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);

                #dSD =nest.GetStatus(self.vnet.spkdetInhibLayers[layer],keys='events')[0];

                # plot spike raster
                ax = plt.subplot(nLayers, 4, (nLayers - layer - 1) * 4 + 3)
                plt.title(layer)
                raster_plot(self.vnet.spkdetInhibLayers[layer])
                # ax.get_yaxis().set_visible(False)
                #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
                #ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);

            plt.show()
