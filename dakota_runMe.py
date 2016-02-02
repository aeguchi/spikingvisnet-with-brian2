from Parameters import *
import sys
import main
import time
import subprocess
import os;
from brian2 import ms

experimentName = "2obj_2"
experimentName = experimentName + time.strftime("_%Y.%m.%d_%H.%M.%S", time.gmtime());
print str(sys.argv);

#set params
outputFile = sys.argv[1];
#taum = float(sys.argv[2])*ms;
taupre = float(sys.argv[2])*ms;
taupost = float(sys.argv[3])*ms;
#conductanceConst_L2L = float(sys.argv[3]);
#dApre = float(sys.argv[3]);

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
#gmax = 0.5;
taum = 20.00*ms;
#ratioPreToPost=5;
#phases = [1,1,1,1,1,2]



#conductanceConst_G2L = float(sys.argv[2]);
#conductanceConst_E2I = float(sys.argv[3]);
#conductanceConst_I2E = float(sys.argv[4]);
#conductanceConst_L2L = float(sys.argv[5]);
#taum = float(sys.argv[6])*ms;


#running simulation
main.loadParams(globals());
main.runSimulation();

#export params
f = open(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/params.txt","w");
f.write(str(sys.argv));
f.close();
source = os.path.split(os.path.realpath(__file__))[0] +"/Dakota/params.in";
destination = os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/params.in"
subprocess.call("cp " + source + " " + destination, shell=True);


#copy the results into dakota folder
source = os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/performance.txt";
destination = outputFile;
subprocess.call("cp " + source + " " + destination, shell=True);