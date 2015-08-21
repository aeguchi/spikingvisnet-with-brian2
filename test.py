from brian2 import *
from Parameters import *

net = Network(collect())
 
layerG = [];
for theta in range(0,len(thetaList)):
    layerG.append(PoissonGroup(layerGDim*layerGDim,numpy.random.rand(layerGDim*layerGDim)*40*Hz ))
net.add(layerG)
#tmp = layerG[0]
#M = SpikeMonitor(layerG[0]);
 
 
monitors = []
for theta in range(0,len(thetaList)):
    monitors.append(SpikeMonitor(layerG[theta]));
 
# a simple run would not include the monitors
net.add(monitors)  # manually add the monitors
 
 
net.run(1*second)
#plot(M.t/ms, M.i, '.')
plot(monitors[0].t/ms, monitors[0].i, '.')
show()

