
from Parameters import *
from brian2 import *


class visnet(object):

    """spiking neural network implementation of visnet"""

    def __init__(self):
        super(visnet, self).__init__()

        self.net = Network(collect())

        self.buildLayers()
        self.buildConnectionsBetweenLayers()
        self.buildSpikeMonitors()

    def buildVN(self):

        return self.net, self.layerG, self.layers, self.inhibLayers, \
            self.spikesG, self.spkdetLayers, self.spkdetInhibLayers

    def buildLayers(self):

        # Gabor input layer

        self.layerG = []
        for theta in range(0, len(thetaList)):
            self.layerG.append(PoissonGroup(
                layerGDim * layerGDim,
                numpy.random.rand(layerGDim * layerGDim) * Hz))

        self.net.add(self.layerG)

        # Excitatory layers

        self.layers = []
        for layer in range(0, nLayers):
            self.layers.append(NeuronGroup(layerDim * layerDim, eqn_membran,
                                           threshold='v>Vt', reset='v = Vr',
                                           refractory=2 * ms))
            self.layers[layer].v = 'Vr'  # + rand() * (Vt - Vr)'
            self.layers[layer].ve = 0 * mV
            self.layers[layer].vi = 0 * mV

        self.net.add(self.layers)

        # Inhibitory self.layers

        self.inhibLayers = []
        for layer in range(0, nLayers):
            self.inhibLayers.append(NeuronGroup(
                layerDim * layerDim, eqn_membran, threshold='v>Vt',
                reset='v = Vr', refractory=2 * ms))
            self.inhibLayers[layer].v = 'Vr'  # + rand() * (Vt - Vr)'
            self.inhibLayers[layer].ve = 0 * mV
            self.inhibLayers[layer].vi = 0 * mV
        self.net.add(self.inhibLayers)

    def buildConnectionsBetweenLayers(self):

        # Connecting neurons in Gabor input layer and neurons in the first
        # ExcitLayer

        connGtoInput = []
        for theta in range(0, len(thetaList)):

            connGtoInput.append(
                Synapses(self.layerG[theta], self.layers[0], pre='ve += 3*we'))
            connGtoInput[theta].connect(True, p=0.2)

        self.net.add(connGtoInput)

        # Connecting neurons between excitatory layers

        connBottomUp = []
        connTopDown = []
        for layer in range(0, nLayers - 1):
            connBottomUp.append(Synapses(
                self.layers[layer], self.layers[layer + 1],
                eqs_stdpSyn, eqs_stdpPre, eqs_stdpPost))

            connBottomUp[layer].connect(True, p=0.1)
            connBottomUp[layer].w[:, :] = rand() * Apre
            connBottomUp[layer].delay[:, :] = rand() * 10 * ms

            connTopDown.append(Synapses(
                self.layers[layer + 1], self.layers[layer],
                eqs_stdpSyn, eqs_stdpPre, eqs_stdpPost))

            connTopDown[layer].connect(True, p=0.1)
            connTopDown[layer].w[:, :] = rand() * Apre
            connTopDown[layer].delay[:, :] = rand() * 10 * ms

        self.net.add(connBottomUp)
        self.net.add(connTopDown)

        # Connecting neurons within layers (excitatory to inhibitory and vv)

        connExIn = []
        connInEx = []
        connRecIn = []
        connRecEx = []

        for layer in range(0, nLayers):
            connExIn.append(
                Synapses(self.layers[layer], self.inhibLayers[layer],
                         'w:1', pre='ve += we'))

            connExIn[layer].connect(True, p=0.2)
            connExIn[layer].delay[:, :] = rand() * 2 * ms

            connInEx.append(
                Synapses(self.inhibLayers[layer], self.layers[layer],
                         'w:1', pre='vi += 3*wi'))

            connInEx[layer].connect(True, p=0.2)
            connInEx[layer].delay[:, :] = rand() * 2 * ms
        self.net.add(connExIn)
        self.net.add(connInEx)

    def buildSpikeMonitors(self):

        self.spikesG = []
        for theta in range(0, len(thetaList)):
            self.spikesG.append(SpikeMonitor(self.layerG[theta]))
        self.net.add(self.spikesG)

        # tmp = self.layerG[0];
        # testSpikes = SpikeMonitor(tmp);

        self.spkdetLayers = []
        for layer in range(0, nLayers):
            self.spkdetLayers.append(SpikeMonitor(self.layers[layer]))
        self.net.add(self.spkdetLayers)

        self.spkdetInhibLayers = []
        for layer in range(0, nLayers):
            self.spkdetInhibLayers.append(
                SpikeMonitor(self.inhibLayers[layer]))
        self.net.add(self.spkdetInhibLayers)
