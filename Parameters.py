# date: 24/11/15
# author: Akihiro Eguchi
# description: a class to define parameters used during simulation

import numpy as np
from cmath import sqrt, log
from brian2 import ms, mV, Hz

#class Parameters:
# nStim = 2
# nTrans = 8

mpdePlotShow = 0; #0: continuous, 1: interrupt
plotGaborAtTraining = False;
plotActivitiesAtTraining = False;
plotWeightsAtTraining = True;

#experimentName = "BO_single"
#imageFolder = "BO_single_imgs"
experimentName = "2obj3"
imageFolder = "simpleImages2obj"

#STRUCTURE OF THE NETWORK
nLayers = 2
layerGDim = 20;
layerDim = 20   # size of layer1 in neurons
inhibLayerDim = 10;

#simulationTime = 100;
trainingTime = 10# * ms;
trainingEpochs = 10;
testingTime = 1000# * ms;

#Param for Filtering:
#ksize = 31#31  # 31 the size of the Gabor kernel. If ksize = (a, b), we then have a Gabor kernel of size a x b pixels. As with many other convolution kernels, ksize is preferably odd and the kernel is a square (just for the sake of uniformity).
bw = 1.5 #spatial bandwidth:
thetaList = [0, np.pi/4, np.pi/2, np.pi*3/4]  #0 #orientation: the orientation of the normal to the parallel stripes of the Gabor function.
lamdaList = [2.0]; #wavelength: the wavelength of the sinusoidal factor in the above equation.
gamma = 0.5 #aspect ratio:  the spatial aspect ratio.
#psiList = [0, np.pi, -np.pi/2, np.pi/2]#phase shift:  phase offset
psiList = [0]#phase shift:  phase offset
sigma = 0.5 # the standard deviation of the Gaussian function used in the Gabor filter.
paddingColor = 128;



#neuron params
taum = 10*ms
taue = 5*ms;
taui = 10*ms;
Vt = -54*mV
Vr = -60*mV
El = -60*mV#-49*mV
refractoryPeriod = 2*ms;
Rmax = 30 * Hz  # max FR
Ee = 0*mV;


# eqn_membran = '''
# dv/dt  = (ve+vi-(v-El))/taum : volt (unless refractory)
# dve/dt = -ve/taue : volt (unless refractory) #incoming excitatory voltage
# dvi/dt = -vi/taui : volt (unless refractory) #incoming inhibitory voltage
# '''

eqn_membran = '''
dv/dt = ((ge + gi) * (Ee-Vr) + El - v) / taum : volt (unless refractory)
dge/dt = -ge/taue : 1 #incoming excitatory voltage
dgi/dt = -gi/taui : 1 #incoming inhibitory voltage
'''



#synapse params:  from < 100 microseconds in very short axons to > 100 ms in very long non-myelinated central axons. 
delayRandOn = True;
weightNormalizationOn = True;
typeOfWeightNormalization = 2; #1:normal, 2:unit norm.(feature rescaling)
delayConst_G2Input = 20*ms;
delayConst_connBottomUp = 20*ms;
delayConst_connExIn = 20*ms;
delayConst_connInEx = 20*ms;
delayConst_connExBind = 20*ms;

#nConnections_connGtoInput = 50;
nConnections_connGtoInput = 10;
fanInRadSigma_connGtoInput = 1;
nConnections_connBottomUp = 50;
nConnections_connExIn = 10;
nConnections_connInEx = 10;
nConnections_connExBind = 10;

#pConnections_connGtoInput = float(nConnections_connGtoInput)/(layerGDim*layerGDim);
pConnections_connBottomUp = float(nConnections_connBottomUp)/(layerDim*layerDim);
pConnections_connExIn = float(nConnections_connExIn)/(layerDim*layerDim);
pConnections_connInEx = float(nConnections_connInEx)/(inhibLayerDim*inhibLayerDim);
pConnections_connExBind = float(nConnections_connExBind)/(layerDim*layerDim);

#Synaptic Connections from Gabor to Layer
we = (0.27/10) # excitatory synaptic weight
wi = (-4.5/30) # inhibitory synaptic weight

conductanceConst_G2L = 10;
eqs_G2LSyn = '''plastic: boolean (shared)'''
eqs_G2LPre ='''ge += conductanceConst_G2L*we'''

conductanceConst_E2I = 2;
eqs_E2ISyn = '''w:1''';
eqs_E2IPre = '''ge += conductanceConst_E2I*we''';

# conductanceConst_E2E = 1;
# eqs_E2ESyn = '''w:1'''
# eqs_E2EPre = '''ge += conductanceConst_E2E*we''';

conductanceConst_I2E = 1;
eqs_I2ESyn = '''w:1'''
eqs_I2EPre = '''gi += conductanceConst_I2E*wi''';


    
#Synaptic Connections with STDP ; for usage, see http://brian2.readthedocs.org/en/2.0b4/examples/synapses.STDP.html
conductanceConst_L2L = 1;
tau_syn_const = 1;
taupre = 20*ms  * tau_syn_const;
taupost = 20*ms  * tau_syn_const;
# wmax = 20 *mV #200 *mV
#Apre = 3.0 *mV #20*mV
#Apost = -Apre*taupre/taupost*1.05
lRate = 1#0.001
gmax = .05#.01
dApre = 0.1
ratioPreToPost = 1.2;
dApost = -dApre * taupre / taupost * ratioPreToPost
dApost *= gmax
dApre *= gmax


eqs_stdpSyn = '''w : 1
                dApre/dt = -Apre / taupre : 1 (event-driven)
                dApost/dt = -Apost / taupost : 1 (event-driven)
                plastic: boolean (shared)'''
                
eqs_stdpPre ='''
            ge += conductanceConst_L2L*w
            Apre += dApre
            w = clip(w+plastic*Apost*lRate, 0, gmax)
             '''            
eqs_stdpPost ='''
             Apost += dApost
             w = clip(w+plastic*Apre*lRate, 0, gmax)
             '''





conductanceConst_L2B = 2;
gmax_bind = .5

eqs_stdpSyn_bind = '''w : 1
                dApre/dt = -Apre / taupre : 1 (event-driven)
                dApost/dt = -Apost / taupost : 1 (event-driven)
                plastic: boolean (shared)'''
                
eqs_stdpPre_bind ='''
            ge += conductanceConst_L2B*w
            Apre += dApre
            w = clip(w+plastic*Apost*lRate, 0, gmax_bind)
             '''            
eqs_stdpPost_bind ='''
             Apost += dApost
             w = clip(w+plastic*Apre*lRate, 0, gmax_bind)
             '''






#condition to introduce topology
fanInRad = layerDim/3;
connCond = 'sqrt((i%layerDim - j%layerDim)**2 + (i/layerDim - j/layerDim)**2) < fanInRad'

randSeed = 1;

