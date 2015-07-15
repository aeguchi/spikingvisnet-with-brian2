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



simulationTime = 500;

#PARAM FOR NEURONS
I_e = 180.0;         # (nA) external current
tau_m = 20.0;      # (ms) membrane time constant
V_th=-55.      # (mV) firing threshold                  
E_L=-70.       # (mV) resting potential

#C_m = 250.0       # (pF) capacitance
#V_reset = -70.0   # (mV) reset potential
#t_ref = 2.0       # (ms) refractory period

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


# Dictionaries (Creation): 
ndict = {"I_e": I_e, "tau_m": tau_m}

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
    "weights": {"uniform":{'min': 0., 'max': 10.}},
    #"delays" : {"linear" :{"c":0.1,"a":0.2}},
    "delays" : {"uniform" :{'min':0.1,'max':0.2}},
    "allow_autapses":False,
    "allow_multapses" :True
}

connForwardDict = {
    "connection_type":"convergent",
    #"synapse_model": "stdp_synapse",
    "kernel":layerSampleRatio,
    #"number_of_connections": int(layer1Dim*layer1Dim*layerSampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 20.}},
    #"delays" : {"linear" :{"c":0.1,"a":0.2}},
    "delays" : {"uniform" :{'min':0.1,'max':50.0}},
    "allow_autapses":False,
    "allow_multapses" :True}


connBackwardDict = {   "connection_type":"convergent",
    #"synapse_model": "stdp_synapse",
    "kernel":layerSampleRatio,
    #"number_of_connections": int(layer2Dim*layer2Dim*layerSampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 10.}},
    #"delays" : {"linear" :{"c":0.1,"a":0.2}},
    "delays" : {"uniform" :{'min':0.1,'max':50.0}},
    "allow_autapses":False,
    "allow_multapses" :True
}

connExInhibDict = {   "connection_type":"convergent",
    #"kernel":layerSampleRatio,
    "number_of_connections": layer2Dim*layer2Dim,
    "weights": {"uniform":{'min': 0., 'max': 1.}},
    "delays" : {"linear" :{"c":0.1,"a":0.2}},
    "allow_autapses":False,
    "allow_multapses" :False
}

connInhibExDict = {   "connection_type":"convergent",
    #"kernel":layerSampleRatio,
    "number_of_connections": inhibLayer1Dim*inhibLayer1Dim,
    "weights": {"uniform":{'min': -10., 'max': 0.}},
    "delays" : {"linear" :{"c":0.1,"a":0.2}},
    "allow_autapses":False,
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