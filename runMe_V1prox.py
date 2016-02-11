from Parameters import *
import main
from brian2 import ms
import InfoAnalysis
import pylab as plt

experimentName="test_V1_norm2"
#imageFolder = "simpleImagesNoTrans"
imageFolder = "simpleImagesNoTrans"
plotGaborAtTraining = False;
plotActivitiesAtTraining = False;
plotWeightsAtTraining = True;
plotPopulationRateOn = False;

#simulation settings
trainingTime = 200;
trainingEpochs = 500;
testingTime = 5000;

typeOfWeightNormalization = 2;
gmax = 2.0;
gmax_bind = 2.0
conductanceConst_E2I = 2.0#0.4;
conductanceConst_I2E = 0.25;
conductanceConst_G2L = 2.0
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


# ia = InfoAnalysis.InfoAnalysis(globals())
# ia.singleCellInfoAnalysis();
# plt.show();