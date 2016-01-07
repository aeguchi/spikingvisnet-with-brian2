from Parameters import *
import sys
import main
import time
import subprocess
import os;
experimentName = experimentName + time.strftime("_%Y.%m.%d_%H.%M.%S", time.gmtime());

print str(sys.argv);

outputFile = sys.argv[1];
conductanceConst_G2L = float(sys.argv[2]);
conductanceConst_E2I = float(sys.argv[3]);
conductanceConst_I2E = float(sys.argv[4]);
conductanceConst_L2L = float(sys.argv[5]);



main.loadParams(globals());
main.runSimulation();

source = os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/performance.txt";
destination = outputFile;

subprocess.call("cp " + source + " " + destination, shell=True);