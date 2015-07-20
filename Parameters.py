import numpy as np
from cmath import sqrt, log

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



simulationTime = 250;

#PARAM FOR NEURONS
I_e = 180.0;         # (nA) external current
tau_m = 20.0;      # (ms) membrane time constant
V_th=-55.      # (mV) firing threshold                  
E_L=-70.       # (mV) resting potential

#C_m = 250.0       # (pF) capacitance
#V_reset = -70.0   # (mV) reset potential
#t_ref = 2.0       # (ms) refractory period

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
tau_minus = 30.0    #defined in post-synaptic neuron (param in postsynaptic neuron)     
stdp_dict = {'tau_plus':15.0,       #Time constant of STDP window, potentiation in ms (tau_minus = 30.0 defined in post-synaptic neuron (param in postsynaptic neuron))
            'lambda':0.005,         #Step size. positive amplitude of the STDP window (A+ by convention).
            'alpha':1.0,            #Asymmetry parameter (scales depressing increments as alpha*lambda)
            'mu_plus':1.0,          #Weight dependence exponent, potentiation
            'mu_minus':1.0,         #Weight dependence exponent, depression
            'Wmax':100.0            #Maximum allowed weight
            }

#PARAMS FOR CONNECTIONS
maxDelay = 50.
minDelay = 1.

#STRUCTURE OF THE NETWORK
nLayers = 2
layerGDim = 20;
layer1Dim = 20   # size of layer1 in neurons
layer2Dim = 20   # size of layer2 in neurons
inhibLayer1Dim = 10;
inhibLayer2Dim = 10;
extentG=[2.,2.] # size of  layer2 in mm
extent1=[2.,2.] # size of  layer1 in mm
extent2=[2.,2.] # size of  layer2 in mm

gaborSampleRatio = 0.1
layerSampleRatio = 0.05
EISampleRatio = 0.3;
IESampleRatio = 1;


# Dictionaries (Creation): 
ndict = {"I_e": I_e, "tau_m": tau_m, "tau_minus": tau_minus}

layerGDict = {"extent" : extent1, # the size of the layer in mm
    "rows" : layerGDim, # the number of rows in this layer ...
    "columns" : layerGDim, # ... and the number of columns
    "elements" : "iaf_neuron"   # the element at each (x,y) coordinate in the grid
    }

layersDict = []
#layer 1
layersDict.append({"extent" : extent1, # the size of the layer in mm
    "rows" : layer1Dim, # the number of rows in this layer ...
    "columns" : layer1Dim, # ... and the number of columns
    "elements" : "iaf_neuron"}) # the element at each (x,y) coordinate in the grid
#layer 2
layersDict.append({"extent" : extent2, # the size of the layer in mm
    "rows" : layer2Dim, # the number of rows in this layer ...
    "columns" : layer2Dim, # ... and the number of columns
    "elements" : "iaf_neuron"}) # the element at each (x,y) coordinate in the grid

inhibLayersDict = []
#inhibLayer 1
inhibLayersDict.append({"extent" : extent1, # the size of the layer in mm
    "rows" : inhibLayer1Dim, # the number of rows in this layer ...
    "columns" : inhibLayer1Dim, # ... and the number of columns
    "elements" : "iaf_neuron"}) # the element at each (x,y) coordinate in the grid
#inhibLayer 2
inhibLayersDict.append({"extent" : extent2, # the size of the layer in mm
    "rows" : inhibLayer2Dim, # the number of rows in this layer ...
    "columns" : inhibLayer2Dim, # ... and the number of columns
    "elements" : "iaf_neuron"}) # the element at each (x,y) coordinate in the grid



# Dictionaries (COnnection):
connGDict = {   "connection_type":"convergent",
    #"synapse_model": "stdp_synapse",
    "kernel":gaborSampleRatio,
    #"number_of_connections": int(layerGDim*layerGDim*gaborSampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 50.}},
    #"delays" : {"linear" :{"c":0.1,"a":0.2}},
    "delays" : {"uniform" :{'min':minDelay,'max':maxDelay}},
    "allow_autapses":False,
    "allow_multapses" :True
}

connForwardDict = {
    "connection_type":"convergent",
    "synapse_model": "stdp_synapse",
    "kernel":layerSampleRatio,
    #"number_of_connections": int(layer1Dim*layer1Dim*layerSampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 100.}},
    #"delays" : {"linear" :{"c":0.1,"a":0.2}},
    "delays" : {"uniform":{'min': minDelay, 'max': maxDelay}},#2.0,#{"uniform" :{'min':0.1,'max':50.0}},
    "allow_autapses":False,
    "allow_multapses" :True}


connBackwardDict = {   "connection_type":"convergent",
    "synapse_model": "stdp_synapse",
    "kernel":layerSampleRatio,
    #"number_of_connections": int(layer2Dim*layer2Dim*layerSampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 50.}},
    #"delays" : {"linear" :{"c":0.1,"a":0.2}},
    "delays" : {"uniform":{'min': minDelay, 'max': maxDelay}},#2.0,# {"uniform" :{'min':0.1,'max':50.0}},
    "allow_autapses":False,
    "allow_multapses" :True
}

connExInhibDict = {   "connection_type":"convergent",
    #"kernel":layerSampleRatio,
    "number_of_connections": int(layer2Dim*layer2Dim*EISampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 10.}},
    "delays" : {"uniform":{'min': minDelay, 'max': maxDelay}},#2.0,#{"linear" :{"c":0.1,"a":0.2}},
    "allow_autapses":False,
    "allow_multapses" :False
}

connInhibExDict = {   "connection_type":"convergent",
    #"kernel":layerSampleRatio,
    "number_of_connections": int(inhibLayer1Dim*inhibLayer1Dim*IESampleRatio),
    "weights": {"uniform":{'min': -8., 'max': 0.}},
    "delays" : {"uniform":{'min': minDelay, 'max': maxDelay}},#2.0,#{"linear" :{"c":0.1,"a":0.2}},
    "allow_autapses":False,
    "allow_multapses" :False
}

connInhibRecDict = {   "connection_type":"convergent",
    #"kernel":layerSampleRatio,
    "number_of_connections": inhibLayer1Dim*inhibLayer1Dim,
    "weights": {"uniform":{'min': 0., 'max': 0.1}},
    "delays" : {"uniform":{'min': minDelay, 'max': 5.}},#5.0,#{"linear" :{"c":0.1,"a":0.2}},
    "allow_autapses":True,
    "allow_multapses" :False
}


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
# }