# date: 24/11/15
# author: Akihiro Eguchi
# description: a class to define parameters used during simulation

import numpy as np
from cmath import sqrt, log
from brian2 import ms, mV, Hz

#class Parameters:
# nStim = 2
# nTrans = 8

#experimentName = "BO_single"
#imageFolder = "BO_single_imgs"
experimentName = "simpleImages2"
imageFolder = "simpleImages"

#STRUCTURE OF THE NETWORK
nLayers = 2
layerGDim = 20;
layerDim = 20   # size of layer1 in neurons
inhibLayerDim = 10;

#simulationTime = 100;
trainingTime = 100# * ms;
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
taum = 20*ms
taue = 5*ms;
taui = 10*ms;
Vt = -50*mV
Vr = -60*mV
El = -60*mV#-49*mV
refractoryPeriod = 2*ms;
Rmax = 20 * Hz  # max FR

eqn_membran = '''
dv/dt  = (ve+vi-(v-El))/taum : volt (unless refractory)
dve/dt = -ve/taue : volt (unless refractory) #incoming excitatory voltage
dvi/dt = -vi/taui : volt (unless refractory) #incoming inhibitory voltage
'''



#synapse params:  from < 100 microseconds in very short axons to > 100 ms in very long non-myelinated central axons. 
delayRandOn = False;
delayConst_connBottomUp = 1*ms;
delayConst_connExIn = 1*ms;
delayConst_connInEx = 5*ms;
delayConst_connExBind = 5*ms;

nConnections_connGtoInput = 50;
nConnections_connBottomUp = 50;
nConnections_connExIn = 10;
nConnections_connInEx = 10;
nConnections_connExBind = 5;

pConnections_connGtoInput = float(nConnections_connGtoInput)/(layerGDim*layerGDim);
pConnections_connBottomUp = float(nConnections_connBottomUp)/(layerDim*layerDim);
pConnections_connExIn = float(nConnections_connExIn)/(layerDim*layerDim);
pConnections_connInEx = float(nConnections_connInEx)/(layerDim*layerDim);
pConnections_connExBind = float(nConnections_connExBind)/(layerDim*layerDim);

#Synaptic Connections from Gabor to Layer
we = (60*0.27/10)*mV # excitatory synaptic weight (voltage)
wi = (-20*4.5/10)*mV # inhibitory synaptic weight

conductanceConst_G2L = 7;
eqs_G2LSyn = '''plastic: boolean (shared)'''
eqs_G2LPre ='''ve += conductanceConst_G2L*we'''

conductanceConst_E2I = 3;
eqs_E2ISyn = '''w:1''';
eqs_E2IPre = '''ve += conductanceConst_E2I*we''';

conductanceConst_I2E = 3;
eqs_I2ESyn = '''w:1'''
eqs_I2EPre = '''vi += conductanceConst_I2E*wi''';


    
#Synaptic Connections with STDP ; for usage, see http://brian2.readthedocs.org/en/latest/resources/tutorials/2-intro-to-brian-synapses.html
const = 5;
taupre = 5*ms  * const;
taupost = 20*ms  * const;
wmax = 20 *mV #200 *mV
Apre = 3.0 *mV #20*mV
Apost = -Apre*taupre/taupost*1.05
lRate = 0.1#0.001
conductanceConst_L2L = 1;


eqs_stdpSyn ='''
             w : volt
             dapre/dt = -apre/taupre : volt (event-driven)
             dapost/dt = -apost/taupost : volt (event-driven)
             plastic: boolean (shared)
             '''
eqs_stdpPre ='''
             v_post += conductanceConst_L2L*w
             apre += Apre
             w = clip(w+plastic*apost*lRate, 0, wmax)
             '''
eqs_stdpPost ='''
             apost += Apost
             w = clip(w+plastic*apre*lRate, 0, wmax)
             '''

#condition to introduce topology
fanInRad = layerDim/3;
connCond = 'sqrt((i%layerDim - j%layerDim)**2 + (i/layerDim - j/layerDim)**2) < fanInRad'


