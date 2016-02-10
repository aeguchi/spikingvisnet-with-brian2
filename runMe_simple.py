from Parameters import *
import main
from brian2 import ms

experimentName="test_norm1_itr300_gmax0.25"
#imageFolder = "simpleImagesNoTrans"
imageFolder = "simpleImagesNoTrans"
plotGaborAtTraining = False;
plotActivitiesAtTraining = False;
plotWeightsAtTraining = True;

#simulation settings
trainingTime = 200;
trainingEpochs = 300;
testingTime = 5000;

typeOfWeightNormalization = 1;
gmax = 0.25;
gmax_bind = 0.25
#layer settings
#layerGDim = 20;
#inhibLayerDim = 5;
#layerDim = 10;

#neuron settings
#taum = 20.00*ms;

#Synaptic settings
# lRate = 1;
# gmax = 0.5;
# taupre = 20*ms
# taupost = 20*ms;
# gmax_bind = 0.1;

#connection setitngs
# nConnections_connGtoInput = 10;
# fanInRadSigma_connGtoInput = 1;
# pConnections_connBottomUp = 0.1;
# pConnections_connExIn = 0.1;
# pConnections_connInEx = 0.1;
# pConnections_connExBind = 0.1;




main.loadParams(globals());
main.runSimulation();