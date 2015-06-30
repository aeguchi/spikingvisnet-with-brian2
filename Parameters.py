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
layerDim = 10;


I_e = 180.0;         # (nA) external current
tau_m = 20.0;      # (ms) membrane time constant
V_th=-55.      # (mV) firing threshold                  
E_L=-70.       # (mV) resting potential

#C_m = 250.0       # (pF) capacitance
#V_reset = -70.0   # (mV) reset potential
#t_ref = 2.0       # (ms) refractory period