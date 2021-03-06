from Parameters import *
import sys
import main
import time
import subprocess
import os;
import spikeAnalysis
from brian2 import ms

objectiveFunc = 1; #1:inforAnalysis, 2:PI, 3:both

experimentName="dakota_visnet_BO_imgs5revised_5connections_smallGmax_trans"
imageFolder = "visnet_BO_imgs5revised_train_mod_two";#"visnet_BO_imgs5revised_train_mod_two"
experimentName = experimentName + time.strftime("_%Y.%m.%d_%H.%M.%S", time.gmtime());
print str(sys.argv);

#set params
outputFile = sys.argv[1];
gmax =  float(sys.argv[2]);
tau_syn_const = float(sys.argv[3]);
#conductanceConst_I2E = float(sys.argv[3]);
#simplifyMultiplitude = 10;
conductanceConst_I2Es[0] = float(sys.argv[4]);
conductanceConst_I2Es[1] = float(sys.argv[5]);
conductanceConst_I2Es[2] = float(sys.argv[6]);



trainingTime = 100.0;
trainingEpochs = 10;#int(trainingEpochs/simplifyMultiplitude);
testingTime = 5000.0;

nSynaptiContact = 1;
nLayers = 4;

#lRate = lRate*3;
# fanInRadSigma_connGtoInput = 0.5
# trainingTime = 200;
# dApre = .1
# 
# 
# conductanceConst_G2L = 20;
# nConnections_connBottomUp = 20;
# nConnections_connExBind = 10;
# conductanceConst_E2I = 5;
# conductanceConst_I2E = 5;
# gmax = 0.5;
# gmax_bind = 0.3;
# #gmax = 0.5;
# taum = 20.00*ms;
#ratioPreToPost=5;
#phases = [1,1,1,1,1,2]

plotGaborAtTraining = False;
plotActivitiesAtTraining = False;
plotWeightsAtTraining = False;
plotPopulationRateOn = False;
ReccurentOn = True;



#conductanceConst_G2L = float(sys.argv[2]);
#conductanceConst_E2I = float(sys.argv[3]);
#conductanceConst_I2E = float(sys.argv[4]);
#conductanceConst_L2L = float(sys.argv[5]);
#taum = float(sys.argv[6])*ms;


#running simulation
main.loadParams(globals());
main.runSimulation();

#spikeAnalysis.loadParams(globals());
#spikeAnalysis.runSpikeAnalysis(2,2,nLayers,PIcalcOn=True,polyAnalysisOn = False,polyHist = False);

#export params
f = open(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/params.txt","w");
f.write(str(sys.argv));
f.close();
source = os.path.split(os.path.realpath(__file__))[0] +"/Dakota/params.in";
destination = os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/params.in"
subprocess.call("cp " + source + " " + destination, shell=True);

#copy the results into dakota folder
if objectiveFunc==1:
    source = os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/performance.txt";
elif objectiveFunc==2:
    source = os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/PI_improvement.txt";
elif objectiveFunc==3:
    source = os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/infoAndPI.txt";
    
destination = outputFile;
subprocess.call("cp " + source + " " + destination, shell=True);