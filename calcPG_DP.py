from brian2 import ms;
import numpy as np;
import pylab as plt;
import pickle
import os;
import errno
import InfoAnalysis
import Queue;
import copy;
import itertools

def loadParams(borrowed_globals):
    globals().update(borrowed_globals);

def runCalcPG(analysisLayer):

    fig = plt.figure();
    #init
    SpikeDataFileName = [str(trainingEpochs+1)+'_spikes_e.pkl', str(trainingEpochs+1)+'_spikes_i.pkl', str(trainingEpochs+1)+'_pikes_b.pkl', str(trainingEpochs+1)+'_pikes_g.pkl'];
    
    f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/" + SpikeDataFileName[0], "rb");
    spikes_e = pickle.load(f);
    f.close();
    #netState_L2L = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/"+ str(trainingEpochs+1) + "_netStates_L2L.pkl", "rb"));
    
    f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/0_netStates_E2E.pkl", "rb");
    netState_E2E_blank = pickle.load(f);
    f.close();
    
    
    f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/"+ str(trainingEpochs+1) + "_netStates_E2E.pkl", "rb");
    netState_E2E = pickle.load(f);
    f.close();
    
    
    drawArrow = True;
    
    
    #analysisLayer = nLayers-1;#(3)
    analysisLayer = analysisLayer-1;
    
    preSynConn = netState_E2E[analysisLayer-1]['i_pre'];
    postSynConn = netState_E2E[analysisLayer-1]['i_post'];
    weightList = netState_E2E[analysisLayer-1]['w'];
    delayList = (netState_E2E[analysisLayer-1]['delay']/ms).astype(int);
    
    nCells = 100;
    sPicked = 3; 
    NbSpikesNeeded = 2;
    NbSpikesMax = 100;
    MaxTimeSpan = 100;
    jitter = 1;
    reflPeriod = jitter+1;
    
    PGCount = 0;
    PGLen = [];
   
    
    
    t_min = 0#840000#15000;#319000;
    t_max = 20000#860000#320000;
            
    #for s_list in [(0, 9, 20)]:#list(itertools.combinations(range(nCells),sPicked)):
    for s_list in list(itertools.combinations(range(nCells),sPicked)):
        #print s_list
        #look for PGs triggered by this combination
        print s_list
                
        
        triggeringConnections = [];
        for i in range(nCells):
            triggeringConnections.append([]);
        
        
        #count connections comming from triggers, to find common triggers output neurons
        for p in s_list:
            cond_pre = preSynConn == p;
            listConnectedFromTp = np.extract(cond_pre, postSynConn);
            for i in listConnectedFromTp:
                if p not in triggeringConnections[i]:#remove duplicate
                    triggeringConnections[i].append(p);
        
        
        for i in range(nCells): 
            PG = []
            PSPTable = np.zeros((nCells,t_max-t_min));
            PSPRecord = np.zeros((nCells,t_max-t_min));
            tLastSpike = np.zeros(nCells)-reflPeriod;
            
        
            if len(triggeringConnections[i]) >= sPicked:
                cond_post = postSynConn == i;
                initSpikeTimings = np.zeros(nCells);
                maxDelay = 0;
                
                #calculate initial spiking time of p
                for p in s_list:
                    cond_pre = preSynConn == p;
                    delay = np.extract(cond_pre & cond_post, delayList);
                    initSpikeTimings[p] = delay[0];#temp focusing only the first elem
                    if maxDelay<delay[0]:
                        maxDelay = delay[0];
                        
                initSpikeTimings = maxDelay-initSpikeTimings;
                
                #update PSPTable based on the every connection from each p
                for p in s_list: 
                    PSPRecord[p,initSpikeTimings[p]]=1;
                    PG.append([p,initSpikeTimings[p]]);
                    cond_pre = preSynConn == p;
                    cellsConnectedFromP = np.extract(cond_pre, postSynConn);
                    delaysConnectedFromP = np.extract(cond_pre, delayList);
                    for cell_i in range(len(cellsConnectedFromP)):
                        timing = initSpikeTimings[p]+delaysConnectedFromP[cell_i];
                        post_cell = cellsConnectedFromP[cell_i];
                        PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+1;
                
                
                
                for t in range(t_max):
                    for p in range(nCells):
                        if ((sum(PSPTable[p,t-jitter:t+1])>=NbSpikesNeeded) and (t>=tLastSpike[p]+reflPeriod)):
                            PSPRecord[p,t]=1;
                            PSPTable[p,t-jitter:t+1] = 0;#reset PSP
                            tLastSpike[p] = t;
                            PG.append([p,t]);
                            
                            cond_pre = preSynConn == p;
                            cellsConnectedFromP = np.extract(cond_pre, postSynConn);
                            delaysConnectedFromP = np.extract(cond_pre, delayList);
                            for cell_i in range(len(cellsConnectedFromP)):
                                timing = t+delaysConnectedFromP[cell_i];
                                post_cell = cellsConnectedFromP[cell_i];
                                PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+1;
                                
                #print PG;                   
                if  len(PG)>4:
                    plt.clf();
                    if drawArrow:
                        plt.hold(True);
                        ax = plt.axes();
                                
                        for t in range(t_max):
                            for p in range(nCells):
                                if (PSPRecord[p][t]==1):
                                    cond_pre = preSynConn == p;
                                    cellsConnectedFromP = np.extract(cond_pre, postSynConn);
                                    delaysConnectedFromP = np.extract(cond_pre, delayList);
                                    for cell_i in range(len(cellsConnectedFromP)):
                                        timing = t+delaysConnectedFromP[cell_i];
                                        post_cell = cellsConnectedFromP[cell_i];
                                        
                                        if (sum(PSPRecord[post_cell,timing:timing+jitter+1])>=1):
                                            dy = (post_cell-p)*0.99;
                                            dx = (delaysConnectedFromP[cell_i])*0.99
                                            ax.arrow(t, p,dx , dy, fc='k', ec='k', head_width=0.05, head_length=0.1)
                
                
                    print s_list
                    print PG;
                    PGTrans = zip(*PG);
                    plt.plot(PGTrans[1], PGTrans[0], 'o', markerfacecolor = 'w');
                    plt.ylabel('Cell Index');
                    plt.xlabel('time [ms]');
                    plt.title('supported PG '+str(s_list));
    #                 plt.hold(True);
                    txt = "Jitter=" + str(jitter) + "ReflPeriod=" + str(reflPeriod) + "\n" +str(PG);
                    fig.text(.1,.1,txt)
                    svname=os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGS_Jitter="+str(jitter)+"_T="+str(s_list)+"_i="+str(i)+"len="+str(len(PG));
                    fig.savefig(svname+".pdf");
                    PGCount+=1;
                    PGLen.append(len(PG));
                    #plt.show();
                    
                    #plt.ylim([0, layerGDim * layerGDim - 1])
                    #plt.xlim([timeBegin, timeBegin+simulationTime])
                    
                    #tmp = input();
                        
    print PGCount;                    
    print mean(PGLen);
