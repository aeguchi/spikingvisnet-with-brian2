from Parameters import *
import sys
import main
import time
import subprocess
import os;
from brian2 import ms

experimentName = experimentName + time.strftime("_%Y.%m.%d_%H.%M.%S", time.gmtime());
print str(sys.argv);

#set params
outputFile = sys.argv[1];
taum = float(sys.argv[2])*ms;
#conductanceConst_L2L = float(sys.argv[3]);
dApre = float(sys.argv[3]);

lRate = 10;
trainingTime = 200;
gmax = 0.5;
#taum = 63.602;
#const




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