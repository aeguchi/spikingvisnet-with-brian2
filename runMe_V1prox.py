from Parameters import *
import main
from brian2 import ms
import InfoAnalysis
import pylab as plt
import spikeAnalysis

experimentName="dakota_visnet_BO_imgs5revised_5connections_2016.04.20_11.25.23"
imageFolder = "visnet_BO_imgs5revised_train_mod_single"

plotGaborAtTraining = False;
plotActivitiesAtTraining = False;
plotWeightsAtTraining = False;
plotPopulationRateOn = False;
ReccurentOn = True;


#set params
#tau_syn_const = float(sys.argv[2]);
gmax =  2.234;
conductanceConst_I2E = 2.564;
#simplifyMultiplitude = 10;

trainingTime = 100.0;
trainingEpochs = 10;#int(trainingEpochs/simplifyMultiplitude);
testingTime = 5000.0;


nSynaptiContact = 5;
nLayers = 4;
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
#main.runSimulation();

spikeAnalysis.loadParams(globals());
spikeAnalysis.runSpikeAnalysis(2,2,nLayers,PIcalcOn=True, polyAnalysisOn = True,polyHist = False);

# spikeAnalysis.loadParams(globals());
# spikeAnalysis.runSpikeAnalysis(2,2,PIcalcOn=True,polyAnalysisOn = False,polyHist = False)

#ia = InfoAnalysis.InfoAnalysis(globals())
#ia.singleCellInfoAnalysis(weightedAnalysis = 1);
#plt.show();