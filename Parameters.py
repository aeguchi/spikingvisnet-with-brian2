import numpy as np
from cmath import sqrt, log

#class Parameters:
nStim = 2
nTrans = 8

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



I_e = 180.0;         # (nA) external current
tau_m = 20.0;      # (ms) membrane time constant
V_th=-55.      # (mV) firing threshold                  
E_L=-70.       # (mV) resting potential

#C_m = 250.0       # (pF) capacitance
#V_reset = -70.0   # (mV) reset potential
#t_ref = 2.0       # (ms) refractory period

layerGDim = 10;
layer1Dim = 10   # size of layer1 in neurons
layer2Dim = 10   # size of layer2 in neurons
extentG=[2.,2.] # size of  layer2 in mm
extent1=[2.,2.] # size of  layer1 in mm
extent2=[2.,2.] # size of  layer2 in mm

gaborSampleRatio = 0.1
layerSampleRatio = 0.05


# Dictionaries: 


ndict = {"I_e": I_e, "tau_m": tau_m}

layerGDict = {"extent" : extent1, # the size of the layer in mm
    "rows" : layerGDim*len(thetaList), # the number of rows in this layer ...
    "columns" : layerGDim, # ... and the number of columns
    "elements" : "iaf_neuron"   # the element at each (x,y) coordinate in the grid
    }

layerGSpkDict = {"extent" : extent1, # the size of the layer in mm
    "rows" : layerGDim*len(thetaList), # the number of rows in this layer ...
    "columns" : layerGDim, # ... and the number of columns
    "elements" : "spike_detector"   # the element at each (x,y) coordinate in the grid
}

layer1Dict = {"extent" : extent1, # the size of the layer in mm
    "rows" : layer1Dim, # the number of rows in this layer ...
    "columns" : layer1Dim, # ... and the number of columns
    "elements" : "iaf_neuron"} # the element at each (x,y) coordinate in the grid

layer1SpkDict = {"extent" : extent1, # the size of the layer in mm
    "rows" : layer1Dim, # the number of rows in this layer ...
    "columns" : layer1Dim, # ... and the number of columns
    "elements" : "spike_detector"} # the element at each (x,y) coordinate in the grid



layer2Dict = {"extent" : extent2, # the size of the layer in mm
    "rows" : layer2Dim, # the number of rows in this layer ...
    "columns" : layer2Dim, # ... and the number of columns
    "elements" : "iaf_neuron"} # the element at each (x,y) coordinate in the grid


layer2SpkDict = {"extent" : extent2, # the size of the layer in mm
    "rows" : layer2Dim, # the number of rows in this layer ...
    "columns" : layer2Dim, # ... and the number of columns
    "elements" : "spike_detector"} # the element at each (x,y) coordinate in the grid

connGDict = {   "connection_type":"convergent",
    "synapse_model": "stdp_synapse",
    "kernel":gaborSampleRatio,
    #"number_of_connections": math.floor(layerGDim*layerGDim*gaborSampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 100.}},
    "delays" : {"linear" :{"c":0.1,"a":0.2}},
    "allow_autapses":False,
    "allow_multapses" :True
}

conn1Dict = {
    "connection_type":"convergent",
    "synapse_model": "stdp_synapse",
    "kernel":layerSampleRatio,
    #"number_of_connections": math.floor(layer1Dim*layer1Dim*layerSampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 100.}},
    "delays" : {"linear" :{"c":0.1,"a":0.2}},
    "allow_autapses":False,
    "allow_multapses" :True}


conn2Dict = {   "connection_type":"convergent",
    "synapse_model": "stdp_synapse",
    "kernel":layerSampleRatio,
    #"number_of_connections": math.floor(layer2Dim*layer2Dim*layerSampleRatio),
    "weights": {"uniform":{'min': 0., 'max': 100.}},
    "delays" : {"linear" :{"c":0.1,"a":0.2}},
    "allow_autapses":False,
    "allow_multapses" :True
}

connSpkDict={"number_of_connections":1, 'connection_type':'convergent'}


spkdetGDict = {"withgid": True,
    "withtime": True,
    "to_memory" :True,
    "to_file" : True,
    "label" : "spkG"
}

spkdet1Dict = {"withgid": True,
    "withtime": True,
    "to_memory" :True,
    "to_file" : True,
    "label" : "spk1"
}

spkdet2Dict = {"withgid": True,
    "withtime": True,
    "to_memory" :True,
    "to_file" : True,
    "label" : "spk2"
}