from Parameters import *
import sys
import main

conductanceConst_G2L = sys.argv[1];
conductanceConst_E2I = sys.argv[2];
conductanceConst_I2E = sys.argv[3];


main.loadParams(globals());
main.runSimulation();