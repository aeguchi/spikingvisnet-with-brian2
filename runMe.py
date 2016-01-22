from Parameters import *
import main
from brian2 import ms

experimentName="ratioPreToPost5"
lRate = 1;
trainingTime = 200;
dApre = .1
nConnections_connBottomUp = 30;
nConnections_connExBind = 20;
conductanceConst_I2E = 5;
gmax = 0.2;
#gmax = 0.5;
#taum = 5.00*ms;
#ratioPreToPost=5;
#phases = [1,1,1,1,1,2]

main.loadParams(globals());
main.runSimulation();