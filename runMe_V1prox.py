from Parameters import *
import main
from brian2 import ms
import InfoAnalysis
import pylab as plt

experimentName="topp_itr1000"
#inhibLayerDim = 10;
#imageFolder = "simpleImagesNoTrans"
imageFolder = "simpleImagesNoTrans"
plotGaborAtTraining = False;
plotActivitiesAtTraining = False;
plotWeightsAtTraining = False;
plotPopulationRateOn = False;
ReccurentOn = True;

#simulation settings
trainingTime = 200.0;
trainingEpochs = 1000;
testingTime = 10000.0;

typeOfWeightNormalization = 1;
type1NormConst = 15.0;

nConnections_connGtoInput = 20;
fanInRadSigma_connGtoInput = 2;

conductanceConst_G2L = 2.5
conductanceConst_E2I = 1.0#0.4;
conductanceConst_I2E = 1.0;
conductanceConst_E2E = 1.0
gmax = 2.0;
gmax_bind = 1.0
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