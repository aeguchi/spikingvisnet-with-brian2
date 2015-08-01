from brian2 import *
from Parameters import *

layerG = [];
for theta in range(0,len(thetaList)):
    layerG.append(PoissonGroup(layerGDim*layerGDim,numpy.random.rand(layerGDim*layerGDim)*40*Hz ))

#tmp = layerG[0]

M = SpikeMonitor(layerG[0]);

run(1*second)
plot(M.t/ms, M.i, '.')

show()
