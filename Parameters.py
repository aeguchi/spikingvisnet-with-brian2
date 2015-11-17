import numpy as np
from cmath import sqrt, log
from brian2 import ms, mV, Hz

#class Parameters:
# nStim = 2
# nTrans = 8

experimentName = "BO_single"
imageFolder = "BO_single_imgs"

#STRUCTURE OF THE NETWORK
nLayers = 2
layerGDim = 20;
layerDim = 20   # size of layer1 in neurons
inhibLayerDim = 10;


#Param for Filtering:
ksize = 31#31  # 31 the size of the Gabor kernel. If ksize = (a, b), we then have a Gabor kernel of size a x b pixels. As with many other convolution kernels, ksize is preferably odd and the kernel is a square (just for the sake of uniformity).
bw = 1.5 #spatial bandwidth:
thetaList = [0, np.pi/4, np.pi/2, np.pi*3/4]  #0 #orientation: the orientation of the normal to the parallel stripes of the Gabor function.
lamdaList = [2.0]; #wavelength: the wavelength of the sinusoidal factor in the above equation.
gamma = 0.5 #aspect ratio:  the spatial aspect ratio.
#psiList = [0, np.pi, -np.pi/2, np.pi/2]#phase shift:  phase offset
psiList = [0]#phase shift:  phase offset
sigma = 0.5 # the standard deviation of the Gaussian function used in the Gabor filter.
paddingColor = 128;



#synapse params
delayConst_connBottomUp = 5*ms;
delayConst_connExIn = 5*ms;
delayConst_connInEx = 5*ms;
delayConst_connExBind = 5*ms;

nConnections_connGtoInput = 50;
nConnections_connBottomUp = 100;
nConnections_connExIn = 10;
nConnections_connInEx = 10;
nConnections_connExBind = 20;

pConnections_connGtoInput = float(nConnections_connGtoInput)/(layerGDim*layerGDim);
pConnections_connBottomUp = float(nConnections_connBottomUp)/(layerDim*layerDim);
pConnections_connExIn = float(nConnections_connExIn)/(layerDim*layerDim);
pConnections_connInEx = float(nConnections_connInEx)/(layerDim*layerDim);
pConnections_connExBind = float(nConnections_connExBind)/(layerDim*layerDim);
#Brian
# eqs ='''
#     dV/dt = (-v + R * I)/tau_m : volt
#     R :ohm
#     I: amp
#     #v_th : volt  # neuron-specific threshold
#     #v_r : volt  # neuron-specific reset
#     '''

taum = 20*ms
taue = 5*ms;
taui = 10*ms;
Vt = -50*mV
Vr = -60*mV
El = -60*mV#-49*mV

eqn_membran = '''
dv/dt  = (ve+vi-(v-El))/taum : volt (unless refractory)
dve/dt = -ve/taue : volt (unless refractory) #incoming excitatory voltage
dvi/dt = -vi/taui : volt (unless refractory) #incoming inhibitory voltage
'''



we = (60*0.27/10)*mV # excitatory synaptic weight (voltage)
wi = (-20*4.5/10)*mV # inhibitory synaptic weight
    
#eqn for STDP ; for usage, see http://brian2.readthedocs.org/en/latest/resources/tutorials/2-intro-to-brian-synapses.html


const = 10;
taupre = taupost = 20*ms  * const;
wmax = 20 *mV #200 *mV
Apre = 3.0 *mV #20*mV
Apost = -Apre*taupre/taupost*1.05
lRate = 0.01#0.001


eqs_stdpSyn ='''
             w : volt
             dapre/dt = -apre/taupre : volt (event-driven)
             dapost/dt = -apost/taupost : volt (event-driven)
             plastic: boolean (shared)
             '''
eqs_stdpPre ='''
             v_post += w
             apre += Apre
             w = clip(w+plastic*apost*lRate, 0, wmax)
             '''
eqs_stdpPost ='''
             apost += Apost
             w = clip(w+plastic*apre*lRate, 0, wmax)
             '''

connCond = '''
        sqrt((i%layerDim - j%layerDim)**2 + (i/layerDim - j/layerDim)**2) < 0.1 * layerDim
        '''


simulationTime = 250;

#PARAM FOR NEURONS
# I_e = 180.0 * nA;         # (nA) external current
# tau_m = 20.0 * ms;      # (ms) membrane time constant
# V_th=-55. * mV      # (mV) firing threshold
# E_L=-70. * mV     # (mV) resting potential
Rmax = 20 * Hz  # max FR


#PARAMS FOR CONNECTIONS
maxDelay = 50.
minDelay = 1.


