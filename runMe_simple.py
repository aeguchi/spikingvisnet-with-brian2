from Parameters import *
import main
from brian2 import ms

experimentName="poly1"
#imageFolder = "simpleImagesNoTrans"
imageFolder = "simpleImagesNoTrans"
plotActivities = 1;
lRate = 1;
fanInRadSigma_connGtoInput = 0.5
trainingTime = 200;
dApre = .1

conductanceConst_G2L = 20;
nConnections_connBottomUp = 20;
nConnections_connExBind = 10;
conductanceConst_E2I = 5;
conductanceConst_I2E = 5;
gmax = 0.5;
gmax_bind = 0.3;
taupre = 100*ms;#200*ms
taupost = 100*ms;#200*ms
#gmax = 0.5;
taum = 20.00*ms;
#ratioPreToPost=5;
#phases = [1,1,1,1,1,2]

main.loadParams(globals());
main.runSimulation();