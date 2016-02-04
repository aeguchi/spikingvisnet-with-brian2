# date: 24/11/15
# author: Akihiro Eguchi
# description: a class to construct the network structure


#from Parameters import *
import brian2 as br
from brian2 import Hz, mV, ms  # explicitly import for readability
import numpy as np


class visnet(object):

    """spiking neural network implementation of visnet"""

    def __init__(self,borrowed_globals):
        super(visnet, self).__init__()
        globals().update(borrowed_globals);
        np.random.seed(randSeed)

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
                               refractory=refractoryPeriod  # length of refractory period
                               )
            )

            self.layers[layer].v = 'Vr'  # + br.rand() * (Vt - Vr)'
            self.layers[layer].ge = 0
            self.layers[layer].gi = 0

        self.net.add(self.layers)

        # Inhibitory layers

        self.inhibLayers = []
        for layer in range(0, nLayers):
            self.inhibLayers.append(br.NeuronGroup(
                inhibLayerDim * inhibLayerDim,
                eqn_membran,
                threshold='v>Vt',
                reset='v = Vr',
                refractory=2 * ms
            )
            )

            self.inhibLayers[layer].v = 'Vr'  # + br.rand() * (Vt - Vr)'
            self.inhibLayers[layer].ge = 0
            self.inhibLayers[layer].gi = 0
        self.net.add(self.inhibLayers)
        
        
        #binding layer
        self.bindingLayer = br.NeuronGroup(
            layerDim * layerDim,
            eqn_membran,
            threshold='v>Vt',
            reset='v = Vr',
            refractory=2 * ms
        )
        self.bindingLayer.v = 'Vr'  # + br.rand() * (Vt - Vr)'
        self.bindingLayer.ge = 0
        self.bindingLayer.gi = 0      
        self.net.add(self.bindingLayer) 
        
        
        
    def buildConnectionsBetweenLayers(self):

        # Connecting neurons in Gabor input layer and neurons in the first
        # ExcitLayer

        self.connGtoInput = []
        for theta in range(0, len(thetaList)):
            self.connGtoInput.append(
                br.Synapses(self.layerG[theta],
                            self.layers[0],
                            eqs_G2LSyn,
                            pre=eqs_G2LPre
                            )
            )
            
        #self.connGtoInput[theta].connect(connCond, p=pConnections_connGtoInput)
        for cellIndex in range(layerDim*layerDim):
            theta = np.random.randint(len(thetaList));
            i_row = int(cellIndex/layerDim);
            i_col = cellIndex%layerDim;

            j_rows = np.random.normal(i_row*layerGDim/layerDim, fanInRadSigma_connGtoInput, nConnections_connGtoInput).astype(int)%layerGDim;
            j_cols = np.random.normal(i_col*layerGDim/layerDim, fanInRadSigma_connGtoInput, nConnections_connGtoInput).astype(int)%layerGDim;
            
            #j_row = int(i_row*layerGDim/layerDim);
            #j_col = int(i_col*layerGDim/layerDim);
            GCellsIndex = layerGDim*j_rows + j_cols;
            self.connGtoInput[theta].connect(GCellsIndex,cellIndex);
            #self.connGtoInput[theta].connect('i==cellIndex', p=pConnections_connGtoInput)
            #self.connGtoInput[theta].connect(connCond2, p=pConnections_connGtoInput)

        #self.connGtoInput[theta].connect(True, p=prob_GtoInput)

        #self.connGtoInput[theta].delay[:, :] = br.rand() * delayConst_G2Input if delayRandOn else delayConst_G2Input
        self.connGtoInput[theta].delay[:, :] = br.rand(len(self.connGtoInput[theta].delay[:, :])) * delayConst_G2Input if delayRandOn else delayConst_G2Input
        self.net.add(self.connGtoInput)

        # Connecting neurons between excitatory layers

        self.connBottomUp = []
        #self.connTopDown = []

        for layer in range(nLayers - 1):

            self.connBottomUp.append(br.Synapses(
                self.layers[layer],
                self.layers[layer + 1],
                eqs_stdpSyn,
                pre = eqs_stdpPre,
                post = eqs_stdpPost
            )
            )


            for cellIndex in range(layerDim*layerDim):
                self.connBottomUp[layer].connect('i==cellIndex', p=pConnections_connBottomUp)
            #self.connBottomUp[layer].connect(True, p=prob_connBottomUp)
            #self.connBottomUp[layer].w[:, :] = br.rand() * gmax
            self.connBottomUp[layer].w[:, :] = br.rand(len(self.connBottomUp[layer].w[:, :])) * gmax
            #self.connBottomUp[layer].delay[:, :] = br.rand() * delayConst_connBottomUp if delayRandOn else delayConst_connBottomUp
            self.connBottomUp[layer].delay[:, :] = br.rand(len(self.connBottomUp[layer].delay[:, :])) * delayConst_connBottomUp if delayRandOn else delayConst_connBottomUp

#             self.connTopDown.append(br.Synapses(
#                 self.layers[layer + 1],
#                 self.layers[layer],
#                 eqs_stdpSyn,
#                 eqs_stdpPre,
#                 eqs_stdpPost
#             )
#             )
# 
#             self.connTopDown[layer].connect(True, p=0.1)
#             self.connTopDown[layer].w[:, :] = br.rand() * Apre
#             self.connTopDown[layer].delay[:, :] = br.rand() * 10 * ms

        self.net.add(self.connBottomUp)
        # self.net.add(self.connTopDown)

        # Connecting neurons within layers (excitatory to inhibitory and vv)

        self.connExIn = []
        self.connInEx = []
        # connRecIn = []
        #connRecEx = []

        for layer in range(0, nLayers):
            #Excitatory -> Inhibitory
            self.connExIn.append(
                br.Synapses(self.layers[layer], self.inhibLayers[layer],
                            eqs_E2ISyn, pre=eqs_E2IPre))
            for cellIndex in range(inhibLayerDim*inhibLayerDim):
                self.connExIn[layer].connect('j==cellIndex', p=pConnections_connExIn)
            #self.connExIn[layer].connect(True, p=prob_connExIn)
            #self.connExIn[layer].delay[:, :] = br.rand() * delayConst_connExIn if delayRandOn else delayConst_connExIn
            self.connExIn[layer].delay[:, :] = br.rand(len(self.connExIn[layer].delay[:, :])) * delayConst_connExIn if delayRandOn else delayConst_connExIn


            #Inhibitory -> Excitatotry 
            self.connInEx.append(
                br.Synapses(self.inhibLayers[layer], self.layers[layer],
                            eqs_I2ESyn, pre=eqs_I2EPre))

            for cellIndex in range(layerDim*layerDim):
                self.connInEx[layer].connect('j==cellIndex', p=pConnections_connInEx);
            #self.connInEx[layer].connect(True, p=prob_connInEx)
            self.connInEx[layer].delay[:, :] = br.rand(len(self.connInEx[layer].delay[:, :])) * delayConst_connInEx if delayRandOn else delayConst_connInEx
            
            #Excitatory -> Excitatory
#             self.connRecEx.append(
#                 br.Synapses(self.layers[layer], self.layers[layer],
#                             eqs_E2ESyn, pre=eqs_E2EPre))
# 
#             for cellIndex in range(layerDim*layerDim):
#                 self.connRecEx[layer].connect('j==cellIndex', p=pConnections_connRecEx);
#             self.connRecEx[layer].delay[:, :] = br.rand() * delayConst_connRecEx if delayRandOn else delayConst_connRecEx
            
            
        self.net.add(self.connExIn)
        self.net.add(self.connInEx)
#         self.net.add(self.connRecEx)
        
        
        self.connExBind = []
        
        for layer in range(0, nLayers):

            self.connExBind.append(br.Synapses(
                self.layers[layer],
                self.bindingLayer,
                eqs_stdpSyn_bind,
                pre = eqs_stdpPre_bind,
                post = eqs_stdpPost_bind
            )
            )

            for cellIndex in range(layerDim*layerDim):
                self.connExBind[layer].connect('j==cellIndex', p=pConnections_connExBind)
            
            #self.connExBind[layer].connect(True, p=prob_connExBind)
            #self.connExBind[layer].w[:, :] = br.rand() * gmax_bind
            self.connExBind[layer].w[:, :] = br.rand(len(self.connExBind[layer].w[:, :])) * gmax_bind
            #self.connExBind[layer].delay[:, :] = br.rand() * delayConst_connExBind if delayRandOn else delayConst_connExBind
            self.connExBind[layer].delay[:, :] = br.rand(len(self.connExBind[layer].delay[:, :])) * delayConst_connExBind if delayRandOn else delayConst_connExBind
        self.net.add(self.connExBind)


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
        
        self.spkdetBindingLayer = br.SpikeMonitor(self.bindingLayer)
        self.net.add(self.spkdetBindingLayer)

    def setGaborFiringRates(self, res_norm):
        ps = len(psiList);
        ss = len(lamdaList);
        ors = len(thetaList);
        
        for p in range(ps):
            for s in range(ss):
                for o in range(ors):
        #for index_filter in range(0, len(res_norm)):
                    index_filter = p*(ss*ors)+s*ors+o;
                    #tmp = np.array(res_norm[p][s][o]);
                    #print tmp.shape
                    r = np.reshape(np.array(res_norm[p][s][o]),layerGDim * layerGDim);
                    #r = np.reshape(np.mean(tmp, axis=3),(layerGDim * layerGDim))
                    # print r
                    self.layerG[index_filter].rates = r * Rmax
                    # To be fixed
                    # print vnet.layerG[index_filter].rates

    def traceReset(self):
        for layer in range(nLayers):
            self.layers[layer].v = 'Vr'  # + br.rand() * (Vt - Vr)'
            self.layers[layer].ge = 0
            self.layers[layer].gi = 0
            self.inhibLayers[layer].v = 'Vr'  # + br.rand() * (Vt - Vr)'
            self.inhibLayers[layer].ge = 0
            self.inhibLayers[layer].gi = 0
            self.connExBind[layer].Apre[:,:] = 0;
            self.connExBind[layer].Apost[:,:] = 0;
            
        for layer in range(0, nLayers - 1):
            self.connBottomUp[layer].Apre[:,:] = 0;
            self.connBottomUp[layer].Apost[:,:] = 0;

        self.bindingLayer.v = 'Vr'  # + br.rand() * (Vt - Vr)'
        self.bindingLayer.ge = 0
        self.bindingLayer.gi = 0
        print "Trace reset"
        
    def setSynapticPlasticity(self, synapticBool):

#        for connection in [self.connGtoInput, self.connBottomUp,self.connTopDown, self.connExIn, self.connInEx]:
        for connection in [self.connGtoInput, self.connBottomUp,self.connExIn, self.connInEx, self.connExBind]:    
            for syn in self.connBottomUp:
                syn.plastic = synapticBool
                
    def getFiringRateMap(self,lDim,layer,timeBegin,simulationTime):
        FRMap = np.zeros((lDim, lDim))
        spikeTrains = layer.spike_trains();
        for row_tmp in range(lDim):
            for col_tmp in range(lDim):
                index_tmp = row_tmp * lDim + col_tmp
                if (len(spikeTrains[index_tmp]) == 0):
                    condition = spikeTrains[index_tmp] > timeBegin
                else:
                    condition = spikeTrains[index_tmp] > timeBegin * ms
                FRMap[row_tmp][col_tmp] = len(np.extract(condition, spikeTrains[index_tmp]));
        FRMap = FRMap/(float(simulationTime)/1000);
        return FRMap;

    def weightNormalization(self):
        for layer in range(nLayers - 1):
            if typeOfWeightNormalization==1:
                w_tmp = self.connBottomUp[layer].w[:, :];
                w_norm =  np.sqrt(np.sum(np.square(w_tmp)));
                self.connBottomUp[layer].w[:, :] = w_tmp/w_norm*gmax;
            elif typeOfWeightNormalization==2:
                w_tmp = self.connBottomUp[layer].w[:, :];
                self.connBottomUp[layer].w[:, :] = (w_tmp - np.min(w_tmp))/(np.max(w_tmp)-np.min(w_tmp))*gmax;
            
        
        

