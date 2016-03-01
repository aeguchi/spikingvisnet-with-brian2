from Parameters import *
import main
from brian2 import ms
import InfoAnalysis
import pylab as plt
import spikeAnalysis

experimentName="test_fourlayers"
imageFolder = "BO_single"

plotGaborAtTraining = False;
plotActivitiesAtTraining = False;
plotWeightsAtTraining = False;
plotPopulationRateOn = False;


ReccurentOn = True;
topDownOn = False;


inputType=1;
nLayers = 4

#simulation settings
trainingTime = 100.0#200.0;
trainingEpochs = 10#1000;
testingTime = 1000.0#10000.0;

#psiList = [-np.pi/2,np.pi/2]
# typeOfWeightNormalization = 1;
# type1NormConst = 15.0;

# nConnections_connGtoInput = 25#20;
# fanInRadSigma_connGtoInput = 1.0;
# fanInRadSigma_connBottomUp = 4.0;


# conductanceConst_G2L = 2.5
# conductanceConst_E2I = 1.0#0.4;
# conductanceConst_I2E = 1.0;
# conductanceConst_E2E = 1.0
# gmax = 2.5;
# gmax_bind = 1.0
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

spikeAnalysis.loadParams(globals());
spikeAnalysis.runSpikeAnalysis(2,2,PIcalcOn=True,polyAnalysisOn = False,polyHist = False);

# spikeAnalysis.loadParams(globals());
# spikeAnalysis.runSpikeAnalysis(2,2,PIcalcOn=True,polyAnalysisOn = False,polyHist = False)

# ia = InfoAnalysis.InfoAnalysis(globals())
# ia.singleCellInfoAnalysis(phases = ['FR_0_blank.pkl'],numBins=3,weightedAnalysis = 1);
#plt.show();