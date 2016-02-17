# date: 24/11/15
# author: Akihiro Eguchi
# description: a class to define parameters used during simulation

import numpy as np
from cmath import sqrt, log
from brian2 import ms, mV, Hz

#class Parameters:
# nStim = 2
# nTrans = 8

modePlotShow = 0; #0: continuous, 1: interrupt
plotGaborAtTraining = False;
plotActivitiesAtTraining = False;
plotWeightsAtTraining = True;
plotPopulationRateOn = False;

ratioTakenToCalcFR = 0.5;

#experimentName = "BO_single"
#imageFolder = "BO_single_imgs"
experimentName = "2obj3"
imageFolder = "simpleImages2obj"

#STRUCTURE OF THE NETWORK
nLayers = 2
layerGDim = 20;
layerDim = 10   # size of layer1 in neurons
inhibLayerDim = 5;
topDownOn = False;
ReccurentOn = False;

#simulationTime = 100;
trainingTime = 10.0# * ms;
trainingEpochs = 10;
testingTime = 1000.0# * ms;

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
Rmax = 30.0 * Hz  # max FR


# #neuron params
# taum = 20.00*ms;
# taue = 2.0*ms;
# taui = 5.0*ms;
# 
# Vt = -53.0*mV #firing threshold potential
# Vr = -74.0*mV #resting potential
# #El = -60.0*mV #-49*mV
# refractoryPeriod = 2*ms;
# 
# Erev_e = 0.0*mV; #excitatory reversal potential
# Erev_i = -70.0*mV; #inhibitory reversal potential

# eqn_membran = '''
# dv/dt  = (ve+vi-(v-El))/taum : volt (unless refractory)
# dve/dt = -ve/taue : volt (unless refractory) #incoming excitatory voltage
# dvi/dt = -vi/taui : volt (unless refractory) #incoming inhibitory voltage
# '''

#eqn_membran = '''
#dv/dt = ((ge - gi) * (Ee-Vr) + El - v) / taum : volt (unless refractory)
#dge/dt = -ge/taue : 1 #incoming excitatory voltage
#dgi/dt = -gi/taui : 1 #incoming inhibitory voltage
#'''

# eqn_membran = '''
# dv/dt =  ((Vr -v) + (ge) * (Erev_e-Vr)  + (-gi) * (Erev_i-Vr) + Iext) / taum : volt (unless refractory)
# dge/dt = -ge/taue : 1 #incoming excitatory voltage
# dgi/dt = -gi/taui : 1 #incoming inhibitory voltage
# Iext : volt #external current input
# '''

#params for membrane eqn
tau_ex = 2.0*ms;
tau_in = 5.0*ms;
Vrev_ex = 0.0*mV; #excitatory reversal potential
Vrev_in = -70.0*mV; #inhibitory reversal potential
refractoryPeriod = 2*ms;

#Excitatory neuron params
taum_ex = 20.00*ms;
Vth_ex = -53.0*mV #firing threshold potential
V0_ex = -74.0*mV #resting potential

eqn_membranEx = '''
dv/dt =  ((V0_ex -v) + (ge * (Vrev_ex-v)) + (gi * (Vrev_in-v))+ Iext) / taum_ex : volt (unless refractory)
dge/dt = -ge/tau_ex : 1 #incoming excitatory voltage
dgi/dt = -gi/tau_in : 1 #incoming inhibitory voltage
Iext : volt #external current input
'''
#to-do is (Vrev_ex-v) right? was originally Vrev_ex-V0_ex and Vrev_in-V0_ex

#Inhibitory neuron params
taum_in = 12.00*ms;
Vth_in = -53.0*mV #firing threshold potential
V0_in = -82.0*mV #resting potential

eqn_membranIn = '''
dv/dt =  ((V0_in -v) + (ge * (Vrev_ex-v)) + (gi * (Vrev_in-v))+ Iext) / taum_in : volt (unless refractory)
dge/dt = -ge/tau_ex : 1 #incoming excitatory voltage
dgi/dt = -gi/tau_in : 1 #incoming inhibitory voltage
Iext : volt #external current input
'''



#synapse params:  from < 100 microseconds in very short axons to > 100 ms in very long non-myelinated central axons. 
delayRandOn = True;
weightNormalizationOn = True;
typeOfWeightNormalization = 2; #1:normal, 2:unit norm.(feature rescaling)
type1NormConst = 10.0;
delayConst_G2Input = 20.0*ms;
delayConst_connBottomUp = 20.0*ms;
delayConst_connExIn = 20.0*ms;
delayConst_connInEx = 20.0*ms;
delayConst_connExBind = 20.0*ms;

delayConst_connTopDown = 20.0*ms;
delayConst_connRecEx = 20.0*ms;

#nConnections_connGtoInput = 50;
nConnections_connGtoInput = 10;
fanInRadSigma_connGtoInput = 1.0;
nConnections_connBottomUp = 50;
fanInRadSigma_connBottomUp = 2.0;
nConnections_connE2I = 10;
fanInRadSigma_connE2I = 2.0;
nConnections_connI2E = 10;
fanInRadSigma_connI2E = 2.0;
nConnections_connExBind = 10;
fanInRadSigma_connExBind = 2.0;
nConnections_connTopDown = 0;
fanInRadSigma_connTopDown = 2.0;
nConnections_connRecEx = 5;
fanInRadSigma_connRecEx = 2.0;


#pConnections_connBottomUp = 0.1;
#pConnections_connExIn = 0.1;
#inhibRad = layerDim/3;
#pConnections_connInEx = 0.1;
pConnections_connExBind = 0.1;
# pConnections_connBottomUp = float(nConnections_connBottomUp)/(layerDim*layerDim);
# pConnections_connExIn = float(nConnections_connExIn)/(layerDim*layerDim);
# pConnections_connInEx = float(nConnections_connInEx)/(inhibLayerDim*inhibLayerDim);
# pConnections_connExBind = float(nConnections_connExBind)/(layerDim*layerDim);
pConnections_connTopDown = float(nConnections_connTopDown)/(layerDim*layerDim);
#pConnections_connRecEx = float(nConnections_connRecEx)/(layerDim*layerDim);


#Synaptic Connections from Gabor to Layer
eqs_Syn = '''w:1'''
eqs_InPre = '''gi += w''';
eqs_ExPre ='''ge += w'''

#we = (0.27/10) # excitatory synaptic weight
#wi = (4.5/30) # inhibitory synaptic weight
conductanceConst_G2L = 0.4;#20*we;
conductanceConst_E2I = 0.2;#10*we;
conductanceConst_E2E = 0.1;#*we; (rec)
conductanceConst_I2E = 0.075;#0.5*wi;
weightRandOn = False; #for G2L, E2I, E2E, I2E


    
#Synaptic Connections with STDP ; for usage, see http://brian2.readthedocs.org/en/2.0b4/examples/synapses.STDP.html
tau_syn_const = 1.0;
taupre = 20.0*ms  * tau_syn_const;
taupost = 20.0*ms  * tau_syn_const;
# wmax = 20 *mV #200 *mV
#Apre = 3.0 *mV #20*mV
#Apost = -Apre*taupre/taupost*1.05
lRate = 0.1#0.001
gmax = 0.5#.01
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
            ge += w
            Apre += dApre
            w = clip(w+plastic*Apost*lRate, 0, gmax)
             '''            
eqs_stdpPost ='''
             Apost += dApost
             w = clip(w+plastic*Apre*lRate, 0, gmax)
             '''





gmax_bind = 0.1
dApre_bind = 0.1
ratioPreToPost_bind = 1.2;
dApost_bind = -dApre_bind * taupre / taupost * ratioPreToPost_bind
dApost_bind *= gmax_bind
dApre_bind *= gmax_bind

eqs_stdpSyn_bind = '''w : 1
                dApre_bind/dt = -Apre_bind / taupre : 1 (event-driven)
                dApost_bind/dt = -Apost_bind / taupost : 1 (event-driven)
                plastic: boolean (shared)'''
                
eqs_stdpPre_bind ='''
            ge += w
            Apre_bind += dApre_bind
            w = clip(w+plastic*Apost_bind*lRate, 0, gmax_bind)
             '''            
eqs_stdpPost_bind ='''
             Apost_bind += dApost_bind
             w = clip(w+plastic*Apre_bind*lRate, 0, gmax_bind)
             '''






#condition to introduce topology
fanInRad = layerDim/3;
connCond = 'sqrt(0-99(i%layerDim - j%layerDim)**2 + (i/layerDim - j/layerDim)**2) < fanInRad'
#somehow sqrt did not seem to work, instead use one below: (but this is not wrap-around)
connCond = '(((i%layerDim) - (j%layerDim))%layerDim)**2 + (((i/layerDim) - (j/layerDim))%layerDim)**2 < inhibRad**2'
randSeed = 1;

