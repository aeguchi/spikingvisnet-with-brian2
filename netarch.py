
from Parameters import *
import brian2 as br
from brian2 import Hz, mV, ms  # explicitly import for readability
import numpy as np


class visnet(object):

    """spiking neural network implementation of visnet"""

    def __init__(self):
        super(visnet, self).__init__()

        self.net = br.Network(br.collect())

        self.buildLayers()
        self.buildConnectionsBetweenLayers()
        self.buildSpikeMonitors()

    def buildLayers(self):

        # Gabor input layer

        self.layerG = []
        for theta in range(0, len(thetaList)):
            self.layerG.append(br.PoissonGroup(
                layerGDim * layerGDim,
                np.random.rand(layerGDim * layerGDim) * Hz))

        self.net.add(self.layerG)

        # Excitatory layers

        self.layers = []
        for layer in range(0, nLayers):
            self.layers.append(br.NeuronGroup(layerDim * layerDim, eqn_membran,
                                              threshold='v>Vt', reset='v = Vr',
                                              refractory=2 * ms))
            self.layers[layer].v = 'Vr'  # + br.rand() * (Vt - Vr)'
            self.layers[layer].ve = 0 * mV
            self.layers[layer].vi = 0 * mV

        self.net.add(self.layers)

        # Inhibitory layers

        self.inhibLayers = []
        for layer in range(0, nLayers):
            self.inhibLayers.append(br.NeuronGroup(
                layerDim * layerDim, eqn_membran, threshold='v>Vt',
                reset='v = Vr', refractory=2 * ms))
            self.inhibLayers[layer].v = 'Vr'  # + br.rand() * (Vt - Vr)'
            self.inhibLayers[layer].ve = 0 * mV
            self.inhibLayers[layer].vi = 0 * mV
        self.net.add(self.inhibLayers)

    def buildConnectionsBetweenLayers(self):

        # Connecting neurons in Gabor input layer and neurons in the first
        # ExcitLayer

        connGtoInput = []
        for theta in range(0, len(thetaList)):

            connGtoInput.append(
                br.Synapses(self.layerG[theta], self.layers[0],
                            pre='ve += 3*we'))
            connGtoInput[theta].connect(True, p=0.2)

        self.net.add(connGtoInput)

        # Connecting neurons between excitatory layers

        connBottomUp = []
        connTopDown = []
        for layer in range(0, nLayers - 1):
            connBottomUp.append(br.Synapses(
                self.layers[layer], self.layers[layer + 1],
                eqs_stdpSyn, eqs_stdpPre, eqs_stdpPost))

            connBottomUp[layer].connect(True, p=0.1)
            connBottomUp[layer].w[:, :] = br.rand() * Apre
            connBottomUp[layer].delay[:, :] = br.rand() * 10 * ms

            connTopDown.append(br.Synapses(
                self.layers[layer + 1], self.layers[layer],
                eqs_stdpSyn, eqs_stdpPre, eqs_stdpPost))

            connTopDown[layer].connect(True, p=0.1)
            connTopDown[layer].w[:, :] = br.rand() * Apre
            connTopDown[layer].delay[:, :] = br.rand() * 10 * ms

        self.net.add(connBottomUp)
        self.net.add(connTopDown)

        # Connecting neurons within layers (excitatory to inhibitory and vv)

        connExIn = []
        connInEx = []
        connRecIn = []
        connRecEx = []

        for layer in range(0, nLayers):
            connExIn.append(
                br.Synapses(self.layers[layer], self.inhibLayers[layer],
                            'w:1', pre='ve += we'))

            connExIn[layer].connect(True, p=0.2)
            connExIn[layer].delay[:, :] = br.rand() * 2 * ms

            connInEx.append(
                br.Synapses(self.inhibLayers[layer], self.layers[layer],
                            'w:1', pre='vi += 3*wi'))

            connInEx[layer].connect(True, p=0.2)
            connInEx[layer].delay[:, :] = br.rand() * 2 * ms
        self.net.add(connExIn)
        self.net.add(connInEx)

    def buildSpikeMonitors(self):

        # gabor layer

        self.spikesG = []
        for theta in range(0, len(thetaList)):
            self.spikesG.append(br.SpikeMonitor(self.layerG[theta]))
        self.net.add(self.spikesG)

        # tmp = self.layerG[0];
        # testSpikes = SpikeMonitor(tmp);

        # excitatory layers

        self.spkdetLayers = []
        for layer in range(0, nLayers):
            self.spkdetLayers.append(br.SpikeMonitor(self.layers[layer]))
        self.net.add(self.spkdetLayers)

        # inhibitory layers

        self.spkdetInhibLayers = []
        for layer in range(0, nLayers):
            self.spkdetInhibLayers.append(
                br.SpikeMonitor(self.inhibLayers[layer]))
        self.net.add(self.spkdetInhibLayers)

    def setGaborFiringRates(self, res_norm):

        for index_filter in range(0, len(thetaList)):

            r = np.reshape(np.mean(res_norm[index_filter], axis=2),
                           (layerGDim * layerGDim))
            # print r
            self.layerG[index_filter].rates = r * Rmax
            # To be fixed
            # print vnet.layerG[index_filter].rates
