import numpy as np
from cmath import sqrt, log
from brian2 import *

#class Parameters:
nStim = 2
nTrans = 8

#Param for Filtering:
ksize = 31#31  # 31 the size of the Gabor kernel. If ksize = (a, b), we then have a Gabor kernel of size a x b pixels. As with many other convolution kernels, ksize is preferably odd and the kernel is a square (just for the sake of uniformity).
bw = 1.5 #spatial bandwidth:
thetaList = [0, np.pi/4, np.pi/2, np.pi*3/4]  #0 #orientation: the orientation of the normal to the parallel stripes of the Gabor function.
lamdaList = [6.0]; #wavelength: the wavelength of the sinusoidal factor in the above equation.
gamma = 0.5 #aspect ratio:  the spatial aspect ratio.
#psiList = [0, np.pi, -np.pi/2, np.pi/2]#phase shift:  phase offset
psiList = [0]#phase shift:  phase offset
imageFolder = "images"
sigma = 0.5 # the standard deviation of the Gaussian function used in the Gabor filter.


#Brian
# eqs ='''
#     dV/dt = (-v + R * I)/tau_m : volt
#     R :ohm
#     I: amp
#     #v_th : volt  # neuron-specific threshold
#     #v_r : volt  # neuron-specific reset
#     '''
taum = 20*ms
taue = 5*ms
taui = 10*ms
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

taupre = taupost = 20*ms
wmax = 200 *mV
Apre = 20 *mV
Apost = -Apre*taupre/taupost*1.05
lRate = 0.001


eqs_stdpSyn ='''
             w : volt
             dapre/dt = -apre/taupre : volt (event-driven)
             dapost/dt = -apost/taupost : volt (event-driven)
             '''
eqs_stdpPre ='''
             v_post += w
             apre += Apre
             w = clip(w+apost*lRate, 0, wmax)
             '''
eqs_stdpPost ='''
             apost += Apost
             w = clip(w+apre*lRate, 0, wmax)
             '''

connCond = '''
        sqrt((i%layerDim - j%layerDim)**2 + (i/layerDim - j/layerDim)**2) < 0.1 * layerDim
        '''

#'''
#dv/dt = (ge+gi-(v+49*mV))/(20*ms) : volt/s
#dge/dt = -ge/(5*ms) : siemens/s
#dgi/dt = -gi/(10*ms) : siemens/s
#'''
 #need to learn thishttp://brian2.readthedocs.org/en/2.0b4/user/models.html


# eqs = '''
# dv/dt = (I-v)/tau : 1
# I : 1
# tau : second
# '''


simulationTime = 250;

#PARAM FOR NEURONS
# I_e = 180.0 * nA;         # (nA) external current
# tau_m = 20.0 * ms;      # (ms) membrane time constant
# V_th=-55. * mV      # (mV) firing threshold
# E_L=-70. * mV     # (mV) resting potential
Rmax = 20 * Hz  # max FR

#C_m = 250.0 pF       # (pF) capacitance
#V_reset = -70.0 mV  # (mV) reset potential
#t_ref = 2.0 s      # (ms) refractory period

#PARAMS FOR STDP
# Name: stdp_synapse - Synapse type for spike-timing dependent plasticity.
# 
# Description:
# stdp_synapse is a connector to create synapses with spike time dependent plasticity (as defined in [1]).
# Here the weight dependence exponent can be set separately for potentiation and depression.
# 
# Examples:
#     multiplicative STDP [2]  mu_plus = mu_minus = 1.0
#    additive STDP       [3]  mu_plus = mu_minus = 0.0
#    Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
#    van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0 
#    http://synergetics.github.io/nest/stdp__connection_8h_source.html
#
# tau_minus = 30.0 ms   #defined in post-synaptic neuron (param in postsynaptic neuron)
# stdp_dict = {'tau_plus':15.0 ms,       #Time constant of STDP window, potentiation in ms (tau_minus = 30.0 defined in post-synaptic neuron (param in postsynaptic neuron))
#             'lambda':0.005,         #Step size. positive amplitude of the STDP window (A+ by convention).
#             'alpha':1.0,            #Asymmetry parameter (scales depressing increments as alpha*lambda)
#             'mu_plus':1.0,          #Weight dependence exponent, potentiation
#             'mu_minus':1.0,         #Weight dependence exponent, depression
#             'Wmax':100.0            #Maximum allowed weight
#             }

#PARAMS FOR CONNECTIONS
maxDelay = 50.
minDelay = 1.

#STRUCTURE OF THE NETWORK
nLayers = 2
layerGDim = 20;
layerDim = 20   # size of layer1 in neurons
inhibLayerDim = 10;
# extentG=[2.,2.] # size of  layer2 in mm
# extent1=[2.,2.] # size of  layer1 in mm
# extent2=[2.,2.] # size of  layer2 in mm
# 
# gaborSampleRatio = 0.1
# layerSampleRatio = 0.05
# EISampleRatio = 0.3;
# IESampleRatio = 1;



# # Dictionaries (Creation): 
# ndict = {"I_e": I_e, "tau_m": tau_m, "tau_minus": tau_minus}
# 
# layerGDict = {"extent" : extent1, # the size of the layer in mm
#     "rows" : layerGDim, # the number of rows in this layer ...
#     "columns" : layerGDim, # ... and the number of columns
#     "elements" : "iaf_neuron"   # the element at each (x,y) coordinate in the grid
#     }
# 
# layersDict = []
# #layer 1
# layersDict.append({"extent" : extent1, # the size of the layer in mm
#     "rows" : layerDim, # the number of rows in this layer ...
#     "columns" : layerDim, # ... and the number of columns
#     "elements" : "iaf_neuron"}) # the element at each (x,y) coordinate in the grid
# #layer 2
# layersDict.append({"extent" : extent2, # the size of the layer in mm
#     "rows" : layerDim, # the number of rows in this layer ...
#     "columns" : layerDim, # ... and the number of columns
#     "elements" : "iaf_neuron"}) # the element at each (x,y) coordinate in the grid
# 
# inhibLayersDict = []
# #inhibLayer 1
# inhibLayersDict.append({"extent" : extent1, # the size of the layer in mm
#     "rows" : inhibLayerDim, # the number of rows in this layer ...
#     "columns" : inhibLayerDim, # ... and the number of columns
#     "elements" : "iaf_neuron"}) # the element at each (x,y) coordinate in the grid
# #inhibLayer 2
# inhibLayersDict.append({"extent" : extent2, # the size of the layer in mm
#     "rows" : inhibLayerDim, # the number of rows in this layer ...
#     "columns" : inhibLayerDim, # ... and the number of columns
#     "elements" : "iaf_neuron"}) # the element at each (x,y) coordinate in the grid
# 
# 
# 
# # Dictionaries (COnnection):
# connGDict = {   "connection_type":"convergent",
#     #"synapse_model": "stdp_synapse",
#     "kernel":gaborSampleRatio,
#     #"number_of_connections": int(layerGDim*layerGDim*gaborSampleRatio),
#     "weights": {"uniform":{'min': 0., 'max': 50.}},
#     #"delays" : {"linear" :{"c":0.1,"a":0.2}},
#     "delays" : {"uniform" :{'min':minDelay,'max':maxDelay}},
#     "allow_autapses":False,
#     "allow_multapses" :True
# }
# 
# connForwardDict = {
#     "connection_type":"convergent",
#     "synapse_model": "stdp_synapse",
#     "kernel":layerSampleRatio,
#     #"number_of_connections": int(layer1Dim*layer1Dim*layerSampleRatio),
#     "weights": {"uniform":{'min': 0., 'max': 100.}},
#     #"delays" : {"linear" :{"c":0.1,"a":0.2}},
#     "delays" : {"uniform":{'min': minDelay, 'max': maxDelay}},#2.0,#{"uniform" :{'min':0.1,'max':50.0}},
#     "allow_autapses":False,
#     "allow_multapses" :True}
# 
# 
# connBackwardDict = {   "connection_type":"convergent",
#     "synapse_model": "stdp_synapse",
#     "kernel":layerSampleRatio,
#     #"number_of_connections": int(layer2Dim*layer2Dim*layerSampleRatio),
#     "weights": {"uniform":{'min': 0., 'max': 50.}},
#     #"delays" : {"linear" :{"c":0.1,"a":0.2}},
#     "delays" : {"uniform":{'min': minDelay, 'max': maxDelay}},#2.0,# {"uniform" :{'min':0.1,'max':50.0}},
#     "allow_autapses":False,
#     "allow_multapses" :True
# }
# 
# <<<<<<< local
# connSpkDict={"number_of_connections":1, 'connection_type':'convergent'}
# 
# 
# spkdetGDict = {"withgid": True,
#     "withtime": True,
# <<<<<<< local
#     "to_memory" :False,
#     "to_file" : True,
# =======
#     "to_memory" :True,
#     "to_file" : False,
# >>>>>>> other
#     "label" : "spkG"
# =======
# connExInhibDict = {   "connection_type":"convergent",
#     #"kernel":layerSampleRatio,
#     "number_of_connections": int(layerDim*layerDim*EISampleRatio),
#     "weights": {"uniform":{'min': 0., 'max': 10.}},
#     "delays" : {"uniform":{'min': minDelay, 'max': maxDelay}},#2.0,#{"linear" :{"c":0.1,"a":0.2}},
#     "allow_autapses":False,
#     "allow_multapses" :False
# >>>>>>> other
# }
# 
# <<<<<<< local
# spkdet1Dict = {"withgid": True,
#     "withtime": True,
# <<<<<<< local
#     "to_memory" :False,
#     "to_file" : True,
# =======
#     "to_memory" :True,
#     "to_file" : False,
# >>>>>>> other
#     "label" : "spk1"
# =======
# connInhibExDict = {   "connection_type":"convergent",
#     #"kernel":layerSampleRatio,
#     "number_of_connections": int(inhibLayerDim*inhibLayerDim*IESampleRatio),
#     "weights": {"uniform":{'min': -8., 'max': 0.}},
#     "delays" : {"uniform":{'min': minDelay, 'max': maxDelay}},#2.0,#{"linear" :{"c":0.1,"a":0.2}},
#     "allow_autapses":False,
#     "allow_multapses" :False
# >>>>>>> other
# }
# 
# <<<<<<< working copy
# <<<<<<< local
# spkdet2Dict = {"withgid": True,
#     "withtime": True,
# <<<<<<< local
#     "to_memory" :False,
#     "to_file" : True,
# =======
#     "to_memory" :True,
#     "to_file" : False,
# >>>>>>> other
#     "label" : "spk2"
# }=======
# =======
# connInhibRecDict = {   "connection_type":"convergent",
#     #"kernel":layerSampleRatio,
#     "number_of_connections": inhibLayerDim*inhibLayerDim,
#     "weights": {"uniform":{'min': 0., 'max': 0.1}},
#     "delays" : {"uniform":{'min': minDelay, 'max': 5.}},#5.0,#{"linear" :{"c":0.1,"a":0.2}},
#     "allow_autapses":True,
#     "allow_multapses" :False
# }
# 
# >>>>>>> destination

#connSpkDict={"number_of_connections":1, 'connection_type':'convergent'}


# spkdetGDict = {"withgid": True,
#     "withtime": True,
#     "to_memory" :True,
#     "to_file" : False,
#     "label" : "spkG"
# }
# 
# spkdet1Dict = {"withgid": True,
#     "withtime": True,
#     "to_memory" :True,
#     "to_file" : False,
#     "label" : "spk1"
# }
# 
# spkdet2Dict = {"withgid": True,
#     "withtime": True,
#     "to_memory" :True,
#     "to_file" : False,
#     "label" : "spk2"
# }>>>>>>> other
