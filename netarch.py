
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
        self.setSynapticPlasticity(True)

    def buildLayers(self):

        # Gabor input layer

        self.layerG = []
        for theta in range(0, len(thetaList)):
            self.layerG.append(br.PoissonGroup(
                layerGDim * layerGDim,  # N neurons
                np.random.rand(layerGDim * layerGDim) * Hz)  # rates
            )

        self.net.add(self.layerG)

        # Excitatory layers

        self.layers = []
        for layer in range(0, nLayers):
            self.layers.append(
                br.NeuronGroup(layerDim * layerDim,  # N neurons
                               eqn_membran,  # differential equations
                               threshold='v>Vt',  # spike condition
                               reset='v = Vr',  # code to execute on reset
                               refractory=2 * ms  # length of refractory period
                               )
            )

            self.layers[layer].v = 'Vr'  # + br.rand() * (Vt - Vr)'
            self.layers[layer].ve = 0 * mV
            self.layers[layer].vi = 0 * mV

        self.net.add(self.layers)

        # Inhibitory layers

        self.inhibLayers = []
        for layer in range(0, nLayers):
            self.inhibLayers.append(br.NeuronGroup(
                layerDim * layerDim,
                eqn_membran,
                threshold='v>Vt',
                reset='v = Vr',
                refractory=2 * ms
            )
            )

            self.inhibLayers[layer].v = 'Vr'  # + br.rand() * (Vt - Vr)'
            self.inhibLayers[layer].ve = 0 * mV
            self.inhibLayers[layer].vi = 0 * mV
        self.net.add(self.inhibLayers)

    def buildConnectionsBetweenLayers(self):

        # Connecting neurons in Gabor input layer and neurons in the first
        # ExcitLayer

        self.connGtoInput = []
        for theta in range(0, len(thetaList)):

            self.connGtoInput.append(
                br.Synapses(self.layerG[theta],
                            self.layers[0],
                            'plastic: boolean (shared)',
                            pre='ve += 3*we'
                            )
            )

            self.connGtoInput[theta].connect(True, p=0.2)

        self.net.add(self.connGtoInput)

        # Connecting neurons between excitatory layers

        self.connBottomUp = []
        self.connTopDown = []

        for layer in range(0, nLayers - 1):

            self.connBottomUp.append(br.Synapses(
                self.layers[layer],
                self.layers[layer + 1],
                eqs_stdpSyn,
                eqs_stdpPre,
                eqs_stdpPost
            )
            )

            self.connBottomUp[layer].connect(True, p=0.1)
            self.connBottomUp[layer].w[:, :] = br.rand() * Apre
            self.connBottomUp[layer].delay[:, :] = br.rand() * 10 * ms

            self.connTopDown.append(br.Synapses(
                self.layers[layer + 1],
                self.layers[layer],
                eqs_stdpSyn,
                eqs_stdpPre,
                eqs_stdpPost
            )
            )

            self.connTopDown[layer].connect(True, p=0.1)
            self.connTopDown[layer].w[:, :] = br.rand() * Apre
            self.connTopDown[layer].delay[:, :] = br.rand() * 10 * ms

        self.net.add(self.connBottomUp)
        self.net.add(self.connTopDown)

        # Connecting neurons within layers (excitatory to inhibitory and vv)

        self.connExIn = []
        self.connInEx = []
        connRecIn = []
        connRecEx = []

        for layer in range(0, nLayers):
            self.connExIn.append(
                br.Synapses(self.layers[layer], self.inhibLayers[layer],
                            'w:1', pre='ve += we'))

            self.connExIn[layer].connect(True, p=0.2)
            self.connExIn[layer].delay[:, :] = br.rand() * 2 * ms

            self.connInEx.append(
                br.Synapses(self.inhibLayers[layer], self.layers[layer],
                            'w:1', pre='vi += 3*wi'))

            self.connInEx[layer].connect(True, p=0.2)
            self.connInEx[layer].delay[:, :] = br.rand() * 2 * ms
        self.net.add(self.connExIn)
        self.net.add(self.connInEx)

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

    def setSynapticPlasticity(self, synapticBool):

        for connection in [self.connGtoInput, self.connBottomUp,
                           self.connTopDown, self.connExIn, self.connInEx]:

            for syn in self.connBottomUp:
                syn.plastic = synapticBool

