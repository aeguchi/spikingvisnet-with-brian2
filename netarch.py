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
        print " - initializing layers.."
        self.buildLayers()
        print " - initializing synapses..";
        self.buildConnectionsBetweenLayers()
        print " - initializing spike monitors.."
        self.buildSpikeMonitors()
        self.setSynapticPlasticity(True)
    

    def buildLayers(self):

        # Gabor input layer layerG [theta(orientation),psi(phase),lambda(wavelength)]
        self.layerG = []
        ps = len(psiList);
        ss = len(lamdaList);
        ors = len(thetaList);
        for p in range(ps):
            self.layerG.append([]);
            for s in range(ss):
                self.layerG[p].append([])
                for o in range(ors):
                    self.layerG[p][s].append(br.PoissonGroup(
                        layerGDim * layerGDim,  # N neurons
                        np.random.rand(layerGDim * layerGDim) * Hz)  # rates
                    )

        self.net.add(self.layerG);

        # Excitatory layers
        print "  ..excitatory layers.."
        self.layers = []
        for layer in range(0, nLayers):
            self.layers.append(
                br.NeuronGroup(layerDim * layerDim,  # N neurons
                               eqn_membranEx,  # differential equations
                               threshold='v>Vth_ex',  # spike condition
                               reset='''v = V0_ex
                                ge = 0
                                gi = 0''',  # code to execute on reset
                               refractory=refractoryPeriod  # length of refractory period
                               )
            )
            self.layers[layer].v = 'V0_ex'  # + br.rand() * (Vt - Vr)'
            self.layers[layer].ge = 0
            self.layers[layer].gi = 0

        self.net.add(self.layers)

        # Inhibitory layers
        print "  ..inhibitory layers.."
        self.inhibLayers = []
        for layer in range(0, nLayers):
            self.inhibLayers.append(br.NeuronGroup(
                inhibLayerDim * inhibLayerDim,
                eqn_membranIn,  # differential equations
                threshold='v>Vth_in',  # spike condition
                reset='''v = V0_in
                        ge = 0
                        gi = 0''',
                refractory=refractoryPeriod
            )
            )
            self.inhibLayers[layer].v = 'V0_in'  # + br.rand() * (Vt - Vr)'
            self.inhibLayers[layer].ge = 0
            self.inhibLayers[layer].gi = 0
        self.net.add(self.inhibLayers)
        
        
        #binding layer
        print "  ..binding layers.."
        self.bindingLayer = br.NeuronGroup(
            layerDim * layerDim,
            eqn_membranEx,  # differential equations
            threshold='v>Vth_ex',  # spike condition
            reset='''v = V0_ex
                    ge = 0
                    gi = 0''',
            refractory=refractoryPeriod
        )
        self.bindingLayer.v = 'V0_ex'  # + br.rand() * (Vt - Vr)'
        self.bindingLayer.ge = 0
        self.bindingLayer.gi = 0      
        self.net.add(self.bindingLayer) 
        
        
        
    def buildConnectionsBetweenLayers(self):

        # Connecting neurons in Gabor input layer and neurons in the first ExcitLayer
        # Gabor input layer layerG [theta(orientation),psi(phase),lambda(wavelength)]
        print "  ..G2Input Syns.."
        self.connGtoInput = []
        
        ps = len(psiList);
        ss = len(lamdaList);
        ors = len(thetaList);
        for p in range(ps):
            self.connGtoInput.append([]);
            for s in range(ss):
                self.connGtoInput[p].append([])
                for o in range(ors):
                    self.connGtoInput[p][s].append(
                        br.Synapses(self.layerG[p][s][o],
                                    self.layers[0],
                                    eqs_Syn,
                                    pre=eqs_ExPre
                                    )
                    )
            
        for i_post_index in range(layerDim*layerDim):
            o = np.random.randint(ors);
            p = np.random.randint(ps);
            s = np.random.randint(ss);
            
            i_post_row = int(i_post_index/layerDim);
            i_post_col = i_post_index%layerDim;

            i_pre_rows = np.random.normal(i_post_row*layerGDim/layerDim, fanInRadSigma_connGtoInput, nConnections_connGtoInput).astype(int)%layerGDim;
            i_pre_cols = np.random.normal(i_post_col*layerGDim/layerDim, fanInRadSigma_connGtoInput, nConnections_connGtoInput).astype(int)%layerGDim;
            
            i_pre_indexes = layerGDim*i_pre_rows + i_pre_cols;
            self.connGtoInput[p][s][o].connect(i_pre_indexes,i_post_index);

        for p in range(ps):
            for s in range(ss):
                for o in range(ors):
                    self.connGtoInput[p][s][o].delay[:, :] = br.rand(len(self.connGtoInput[p][s][o].delay[:, :])) * delayConst_G2Input if delayRandOn else delayConst_G2Input
                    self.connGtoInput[p][s][o].w[:, :] = br.rand(len(self.connGtoInput[p][s][o].w[:, :])) * conductanceConst_G2L if weightRandOn else conductanceConst_G2L;
        self.net.add(self.connGtoInput)




        # Connecting neurons within layers (excitatory to inhibitory and vv)
        print "  ..ExIn InEx RecEx Syns.."
        self.connExIn = []
        self.connInEx = []
        # connRecIn = []
        if ReccurentOn:
            self.connRecEx = []

        for layer in range(0, nLayers):
            #Excitatory -> Inhibitory
            self.connExIn.append(
                br.Synapses(self.layers[layer], self.inhibLayers[layer],
                            eqs_Syn, pre=eqs_ExPre))
#             for cellIndex in range(inhibLayerDim*inhibLayerDim):
#                 self.connExIn[layer].connect('j==cellIndex', p=pConnections_connExIn)
            for i_post_index in range(inhibLayerDim*inhibLayerDim):
                i_post_row = int(i_post_index/inhibLayerDim);
                i_post_col = i_post_index%inhibLayerDim;
     
                i_pre_rows = np.random.normal(i_post_row*layerDim/inhibLayerDim, fanInRadSigma_connE2I, nConnections_connE2I).astype(int)%layerDim;
                i_pre_cols = np.random.normal(i_post_col*layerDim/inhibLayerDim, fanInRadSigma_connE2I, nConnections_connE2I).astype(int)%layerDim;
                 
                i_pre_indexes = layerDim*i_pre_rows + i_pre_cols;
                self.connExIn[layer].connect(i_pre_indexes,i_post_index);
            
            self.connExIn[layer].w[:, :] = br.rand(len(self.connExIn[layer].w[:, :])) * conductanceConst_E2I if weightRandOn else conductanceConst_E2I
            self.connExIn[layer].delay[:, :] = br.rand(len(self.connExIn[layer].delay[:, :])) * delayConst_connExIn if delayRandOn else delayConst_connExIn


            #Inhibitory -> Excitatotry 
            self.connInEx.append(
                br.Synapses(self.inhibLayers[layer], self.layers[layer],
                            eqs_Syn, pre=eqs_InPre))

            for i_post_index in range(layerDim*layerDim):
                i_post_row = int(i_post_index/layerDim);
                i_post_col = i_post_index%layerDim;
    
                i_pre_rows = np.random.normal(i_post_row*inhibLayerDim/layerDim, fanInRadSigma_connI2E, nConnections_connI2E).astype(int)%inhibLayerDim;
                i_pre_cols = np.random.normal(i_post_col*inhibLayerDim/layerDim, fanInRadSigma_connI2E, nConnections_connI2E).astype(int)%inhibLayerDim;
                
                i_pre_indexes = inhibLayerDim*i_pre_rows + i_pre_cols;
                self.connInEx[layer].connect(i_pre_indexes,i_post_index);
            
            self.connInEx[layer].w[:, :] = br.rand(len(self.connInEx[layer].w[:, :])) * conductanceConst_I2E if weightRandOn else conductanceConst_I2E
            self.connInEx[layer].delay[:, :] = br.rand(len(self.connInEx[layer].delay[:, :])) * delayConst_connInEx if delayRandOn else delayConst_connInEx

            #Excitatory -> Excitatory
            if ReccurentOn:
                self.connRecEx.append(
                    br.Synapses(self.layers[layer], self.layers[layer],
                                eqs_Syn, pre=eqs_ExPre))
     
                for i_post_index in range(layerDim*layerDim):
                    i_row = int(i_post_index/layerDim);
                    i_col = i_post_index%layerDim;
                    i_pre_rows = np.random.normal(i_row, fanInRadSigma_connRecEx, nConnections_connRecEx).astype(int)%layerDim;
                    i_pre_cols = np.random.normal(i_col, fanInRadSigma_connRecEx, nConnections_connRecEx).astype(int)%layerDim;
                    i_pre_indexes = layerDim*i_pre_rows + i_pre_cols;
                    self.connRecEx[layer].connect(i_pre_indexes,i_post_index);
#                 for cellIndex in range(layerDim*layerDim):
#                     self.connRecEx[layer].connect('j==cellIndex', p=pConnections_connRecEx);
                self.connRecEx[layer].w[:, :] = br.rand(len(self.connRecEx[layer].w[:, :])) * conductanceConst_E2E if weightRandOn else conductanceConst_E2E
                self.connRecEx[layer].delay[:, :] = br.rand(len(self.connRecEx[layer].delay[:, :])) * delayConst_connRecEx if delayRandOn else delayConst_connRecEx
    
        self.net.add(self.connExIn)
        self.net.add(self.connInEx)
        if ReccurentOn:
            self.net.add(self.connRecEx)



        
        #STDP conn        
        # Connecting neurons between excitatory layers
        
        print "  ..BottomUp and TopDown Syns.."
        self.connBottomUp = []
        if topDownOn:
            self.connTopDown = []
            
        for layer in range(nLayers - 1):
            self.connBottomUp.append(br.Synapses(
                self.layers[layer],
                self.layers[layer + 1],
                eqs_stdpSyn,
                pre = eqs_stdpPre,
                post = eqs_stdpPost
            )
            )

#             for cellIndex in range(layerDim*layerDim):
#                 self.connBottomUp[layer].connect('i==cellIndex', p=pConnections_connBottomUp)

            for cellIndex in range(layerDim*layerDim):
                i_row = int(cellIndex/layerDim);
                i_col = cellIndex%layerDim;
                j_rows = np.random.normal(i_row, fanInRadSigma_connBottomUp, nConnections_connBottomUp).astype(int)%layerDim;
                j_cols = np.random.normal(i_col, fanInRadSigma_connBottomUp, nConnections_connBottomUp).astype(int)%layerDim;
                PreCellsIndex = layerDim*j_rows + j_cols;
                self.connBottomUp[layer].connect(PreCellsIndex,cellIndex);
            
            self.connBottomUp[layer].w[:, :] = br.rand(len(self.connBottomUp[layer].w[:, :])) * gmax
            self.connBottomUp[layer].delay[:, :] = br.rand(len(self.connBottomUp[layer].delay[:, :])) * delayConst_connBottomUp if delayRandOn else delayConst_connBottomUp
            
            if topDownOn:
                self.connTopDown.append(br.Synapses(
                    self.layers[layer + 1],
                    self.layers[layer],
                    eqs_stdpSyn,
                    pre = eqs_stdpPre,
                    post = eqs_stdpPost
                )
                )
                
#                 for cellIndex in range(layerDim*layerDim):
#                     self.connTopDown[layer].connect('i==cellIndex', p=pConnections_connTopDown)   
                for i_post_index in range(layerDim*layerDim):
                    i_row = int(i_post_index/layerDim);
                    i_col = i_post_index%layerDim;
                    i_pre_rows = np.random.normal(i_row, fanInRadSigma_connTopDown, nConnections_connTopDown).astype(int)%layerDim;
                    i_pre_cols = np.random.normal(i_col, fanInRadSigma_connTopDown, nConnections_connTopDown).astype(int)%layerDim;
                    i_pre_indexes = layerDim*i_pre_rows + i_pre_cols;
                    self.connTopDown[layer].connect(i_pre_indexes,i_post_index);    
                
                self.connTopDown[layer].w[:, :] = br.rand(len(self.connTopDown[layer].w[:, :])) * gmax
                self.connTopDown[layer].delay[:, :] = br.rand(len(self.connTopDown[layer].delay[:, :])) * delayConst_connTopDown if delayRandOn else delayConst_connTopDown
        self.net.add(self.connBottomUp)
        if topDownOn:
            self.net.add(self.connTopDown)
        
        

        #binding layer conn
        print "  ..ExBind Syns.."
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

#             for cellIndex in range(layerDim*layerDim):
#                 self.connExBind[layer].connect('j==cellIndex', p=pConnections_connExBind)

            for i_post_index in range(layerDim*layerDim):
                i_post_row = int(i_post_index/layerDim);
                i_post_col = i_post_index%layerDim;
     
                i_pre_rows = np.random.normal(i_post_row, fanInRadSigma_connExBind, nConnections_connExBind).astype(int)%layerDim;
                i_pre_cols = np.random.normal(i_post_col, fanInRadSigma_connExBind, nConnections_connExBind).astype(int)%layerDim;
                 
                i_pre_indexes = layerDim*i_pre_rows + i_pre_cols;
                self.connExBind[layer].connect(i_pre_indexes,i_post_index);
            
            self.connExBind[layer].w[:, :] = br.rand(len(self.connExBind[layer].w[:, :])) * gmax_bind
            self.connExBind[layer].delay[:, :] = br.rand(len(self.connExBind[layer].delay[:, :])) * delayConst_connExBind if delayRandOn else delayConst_connExBind
        self.net.add(self.connExBind)


    def buildSpikeMonitors(self):

        # gabor layer
        self.spikesG = []
        
        ps = len(psiList);
        ss = len(lamdaList);
        ors = len(thetaList);
        for p in range(ps):
            self.spikesG.append([]);
            for s in range(ss):
                self.spikesG[p].append([])
                for o in range(ors):
                    self.spikesG[p][s].append(br.SpikeMonitor(self.layerG[p][s][o]))
        self.net.add(self.spikesG)


        # excitatory layers
        self.spkdetLayers = []
        if plotPopulationRateOn:
            self.popMonLayers =[]
        for layer in range(0, nLayers):
            self.spkdetLayers.append(br.SpikeMonitor(self.layers[layer]))
            if plotPopulationRateOn:
                self.popMonLayers.append(br.PopulationRateMonitor(self.layers[layer]))
        self.net.add(self.spkdetLayers)
        if plotPopulationRateOn:
            self.net.add(self.popMonLayers)

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
                    #self.layerG[index_filter].rates = r * Rmax
                    self.layerG[p][s][o].rates = r * Rmax
                    
                    # To be fixed
                    # print vnet.layerG[index_filter].rates

    def traceReset(self):
        #init neurons
        for layer in range(nLayers):
            self.layers[layer].v = 'V0_ex'  # + br.rand() * (Vt - Vr)'
            self.layers[layer].ge = 0
            self.layers[layer].gi = 0
            self.inhibLayers[layer].v = 'V0_in'  # + br.rand() * (Vt - Vr)'
            self.inhibLayers[layer].ge = 0
            self.inhibLayers[layer].gi = 0
        
        self.bindingLayer.v = 'V0_ex'  # + br.rand() * (Vt - Vr)'
        self.bindingLayer.ge = 0
        self.bindingLayer.gi = 0
        
        #init syn
        for layer in range(nLayers):
            self.connExBind[layer].Apre_bind[:,:] = 0;
            self.connExBind[layer].Apost_bind[:,:] = 0;
            
        for layer in range(nLayers - 1):
            self.connBottomUp[layer].Apre[:,:] = 0;
            self.connBottomUp[layer].Apost[:,:] = 0;

        print "  ..Trace reset.."
        
    def setSynapticPlasticity(self, synapticBool):

#        for connection in [self.connGtoInput, self.connBottomUp,self.connTopDown, self.connExIn, self.connInEx]:
        for connection in [self.connBottomUp,self.connExBind]:    
            for syn in self.connBottomUp:
                syn.plastic = synapticBool
                
    def getFiringRateMap(self,lDim,layer,timeBegin,simulationTime):
        FRMap = np.zeros((lDim, lDim))
        spikeTrains = layer.spike_trains();
        for row_tmp in range(lDim):
            for col_tmp in range(lDim):
                index_tmp = row_tmp * lDim + col_tmp
                if (len(spikeTrains[index_tmp]) == 0):
                    condition = spikeTrains[index_tmp] > timeBegin + simulationTime*ratioTakenToCalcFR
                else:
                    condition = spikeTrains[index_tmp] > (timeBegin +simulationTime * ratioTakenToCalcFR) * ms
                FRMap[row_tmp][col_tmp] = len(np.extract(condition, spikeTrains[index_tmp]));
        FRMap = FRMap/(float(simulationTime*ratioTakenToCalcFR)/1000);
        return FRMap;

    def weightNormalization(self):
        for layer in range(nLayers - 1):
            if typeOfWeightNormalization==1:
                w_tmp = self.connBottomUp[layer].w[:, :];
                w_norm =  np.sqrt(np.sum(np.square(w_tmp)));
                self.connBottomUp[layer].w[:, :] = w_tmp/w_norm*gmax*type1NormConst;
                
                w_tmp = self.connExBind[layer].w[:, :];
                w_norm =  np.sqrt(np.sum(np.square(w_tmp)));
                self.connExBind[layer].w[:, :] = w_tmp/w_norm*gmax_bind*type1NormConst;
            elif typeOfWeightNormalization==2:
                w_tmp = self.connBottomUp[layer].w[:, :];
                self.connBottomUp[layer].w[:, :] = (w_tmp - np.min(w_tmp))/(np.max(w_tmp)-np.min(w_tmp))*gmax;

                w_tmp = self.connExBind[layer].w[:, :];
                self.connExBind[layer].w[:, :] = (w_tmp - np.min(w_tmp))/(np.max(w_tmp)-np.min(w_tmp))*gmax_bind;
            
        
    def saveStates(self,dest,itr):
        synStates_GtoInput = []
        ps = len(psiList);
        ss = len(lamdaList);
        ors = len(thetaList);
        for p in range(ps):
            synStates_GtoInput.append([]);
            for s in range(ss):
                synStates_GtoInput[p].append([])
                for o in range(ors):
                    synStates_GtoInput[p][s].append(self.connGtoInput[p][s][o].get_states(['i_pre','i_post','delay','w']))
        
        synStates_BottomUp = []
        for layer in range(nLayers-1):
            synStates_BottomUp.append(self.connBottomUp[layer].get_states(['i_pre','i_post','delay','w']));
        
        synStates_InEx = []
        for layer in range(nLayers):
            synStates_InEx.append(self.connInEx[layer].get_states(['i_pre','i_post','delay','w']));
            
        synStates_ExIn = []
        for layer in range(nLayers):
            synStates_ExIn.append(self.connExIn[layer].get_states(['i_pre','i_post','delay','w']));
            
        synStates_ExBind = []
        for layer in range(nLayers):
            synStates_ExBind.append(self.connExBind[layer].get_states(['i_pre','i_post','delay','w']));
        
        if ReccurentOn:
            synStates_ExEx = []
            for layer in range(nLayers):
                synStates_ExEx.append(self.connRecEx[layer].get_states(['i_pre','i_post','delay','w']));
        
        
        pickle.dump(synStates_GtoInput, open(dest + str(itr)+"_netStates_G2L.pkl", "wb"));
        pickle.dump(synStates_BottomUp, open(dest + str(itr)+"_netStates_L2L.pkl", "wb"));
        pickle.dump(synStates_InEx, open(dest + str(itr)+"_netStates_I2E.pkl", "wb"));
        pickle.dump(synStates_ExIn, open(dest + str(itr)+"_netStates_E2I.pkl", "wb"));
        pickle.dump(synStates_ExBind, open(dest + str(itr)+"_netStates_L2B.pkl", "wb"));
        if ReccurentOn:
            pickle.dump(synStates_ExEx, open(dest + str(itr)+"_netStates_E2E.pkl", "wb"));
        
    
    def saveSpikes(self,dest,itr):
        spikes_g = [];
        spikes_e = [];
        spikes_i = [];
        
        
        
        ps = len(psiList);
        ss = len(lamdaList);
        ors = len(thetaList);
        for p in range(ps):
            spikes_g.append([]);
            for s in range(ss):
                spikes_g[p].append([])
                for o in range(ors):
                    tmp_g = self.spikesG[p][s][o].spike_trains();
                    for i in range(layerGDim*layerGDim):
                        if len(tmp_g[i])>0:
                            tmp_g[i] = np.sort(tmp_g[i]/ms);
                    spikes_g[p][s].append(tmp_g);       

            
        for layer in range(nLayers):
            tmp_e = self.spkdetLayers[layer].spike_trains();
            tmp_i = self.spkdetInhibLayers[layer].spike_trains();
            for i in range(layerDim*layerDim):
                if len(tmp_e[i]>0):
                    tmp_e[i] = np.sort(tmp_e[i]/ms);
            for i in range(inhibLayerDim*inhibLayerDim):
                if len(tmp_i[i]>0):
                    tmp_i[i] = np.sort(tmp_i[i]/ms);
            spikes_e.append(tmp_e);
            spikes_i.append(tmp_i);
    
        spikes_b = self.spkdetBindingLayer.spike_trains();
        for i in range(layerDim*layerDim):
            if len(spikes_b[i]>0):
                spikes_b[i] = np.sort(spikes_b[i]/ms);
        
        pickle.dump(spikes_e, open(dest + str(itr)+"_spikes_e.pkl", "wb"))
        pickle.dump(spikes_i, open(dest + str(itr)+"_spikes_i.pkl", "wb"))
        pickle.dump(spikes_b, open(dest + str(itr)+"_spikes_b.pkl", "wb"))
    
            
