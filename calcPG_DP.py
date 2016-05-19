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

    
    #init
    SpikeDataFileName = [str(trainingEpochs+1)+'_spikes_e.pkl', str(trainingEpochs+1)+'_spikes_i.pkl', str(trainingEpochs+1)+'_pikes_b.pkl', str(trainingEpochs+1)+'_pikes_g.pkl'];
    
    f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/" + SpikeDataFileName[0], "rb");
    spikes_e = pickle.load(f);
    f.close();
    #netState_L2L = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/"+ str(trainingEpochs+1) + "_netStates_L2L.pkl", "rb"));
    
    
    f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/"+ str(trainingEpochs+1) + "_netStates_L2L.pkl", "rb");
    netState_L2L = pickle.load(f);
    f.close();
    
    
    f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/0_netStates_E2E.pkl", "rb");
    netState_E2E_blank = pickle.load(f);
    f.close();
    
    
    f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/"+ str(trainingEpochs+1) + "_netStates_E2E.pkl", "rb");
    netState_E2E = pickle.load(f);
    f.close();
    
    
    drawArrow = True;
    plotFigure = False;
    
    nCells = 100;
    sPicked = 3; 
    NbSpikesNeeded = 3;
    NbSpikesMax = 200;
    #MaxTimeSpan = 100;
    jitter = 1;
    reflPeriod = jitter+1;
    spikingTh = 3
    decayRatio = 0.9;
    
    PGType = 1; #1:supported PG, 2:adapted PGs
    
    if plotFigure:
        fig = plt.figure();
    
    #analysisLayer = nLayers-1;#(3)
    analysisLayer = analysisLayer-1;
    
    #net state of one layer before the final output layer
    preSynConn_FF = netState_L2L[analysisLayer-1]['i_pre'];
    postSynConn_FF = netState_L2L[analysisLayer-1]['i_post']+nCells;
    weightList_FF = netState_L2L[analysisLayer-1]['w'];
    delayList_FF = (netState_L2L[analysisLayer-1]['delay']/ms).astype(int);

    preSynConn_RC = netState_E2E[analysisLayer-1]['i_pre']+nCells;
    postSynConn_RC = netState_E2E[analysisLayer-1]['i_post']+nCells;
    weightList_RC = netState_E2E[analysisLayer-1]['w'];
    delayList_RC = (netState_E2E[analysisLayer-1]['delay']/ms).astype(int);
    
    #net state of the final output layer
    preSynConn_RC_o = netState_E2E[analysisLayer-1]['i_pre']+nCells;
    postSynConn_RC_o = netState_E2E[analysisLayer-1]['i_post']+nCells;
    weightList_RC_o = netState_E2E[analysisLayer-1]['w'];
    delayList_RC_o = (netState_E2E[analysisLayer-1]['delay']/ms).astype(int);
    

    
    
    preSynConn = np.concatenate((preSynConn_FF,preSynConn_RC,preSynConn_RC_o));
    postSynConn = np.concatenate((postSynConn_FF,postSynConn_RC,postSynConn_RC_o));
    weightList = np.concatenate((weightList_FF,weightList_RC,weightList_RC_o));
    delayList = np.concatenate((delayList_FF,delayList_RC,delayList_RC_o));
    

    PGCount = 0;
    PGLen = [];
   
    
    
    t_min = 0#840000#15000;#319000;
    t_max = 1000#860000#320000;
            
    #for s_list in [(0, 9, 20)]:#list(itertools.combinations(range(nCells),sPicked)):
    for s_list in list(itertools.combinations(range(nCells),sPicked)):
        #print s_list
        #look for PGs triggered by this combination
        print s_list
                
        
        triggeringConnections = [];
        for i in range(nCells,nCells*2):
            triggeringConnections.append([]);
            
        
        #count connections comming from triggers, to find common triggers output neurons
        for p in s_list:
            cond_pre = preSynConn_FF == p;
            listConnectedFromTp = np.extract(cond_pre, postSynConn_FF);
            for i in listConnectedFromTp:
                if p not in triggeringConnections[i-nCells]:#remove duplicate
                    triggeringConnections[i-nCells].append(p);
        
        
        for i in range(nCells,nCells*2): 
            PG = []
            PSPTable = np.zeros((nCells*2,t_max-t_min));
            PSPRecord = np.zeros((nCells*2,t_max-t_min));
            tLastSpike = np.zeros(nCells*2)-reflPeriod;
            
        
            if len(triggeringConnections[i-nCells]) >= sPicked:
                cond_post = postSynConn_FF == i;
                initSpikeTimings = np.zeros(nCells);
                maxDelay = 0;
                
                #calculate initial spiking time of p
                for p in s_list:
                    cond_pre = preSynConn_FF == p;
                    delay = np.extract(cond_pre & cond_post, delayList_FF);
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
                    if PGType==1:
                        for cell_i in range(len(cellsConnectedFromP)):
                            timing = initSpikeTimings[p]+delaysConnectedFromP[cell_i];
                            post_cell = cellsConnectedFromP[cell_i];
                            PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+1;
                    elif PGType==2:
                        weightsConnectedFromP = np.extract(cond_pre, weightList);
                        for cell_i in range(len(cellsConnectedFromP)):
                            timing = initSpikeTimings[p]+delaysConnectedFromP[cell_i];
                            post_cell = cellsConnectedFromP[cell_i];
                            dPSP = weightsConnectedFromP[cell_i]
                            PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+dPSP;
                 
                for t in range(t_min,t_max):
                    if len(PG)<NbSpikesMax:
                        for p in range(nCells*2):
                            if PGType==1:
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
                                        if timing>=t_max:
                                            continue;
                                        post_cell = cellsConnectedFromP[cell_i];
                                        PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+1;
                            elif PGType==2:
                                if t>0:
                                    PSPTable[p,t]=PSPTable[p,t-1]*decayRatio+PSPTable[p,t];
                                if ((PSPTable[p,t]>=spikingTh) and (t>=tLastSpike[p]+reflPeriod)):
                                    PSPRecord[p,t]=1;
                                    PSPTable[p,t] = 0;#reset PSP
                                    tLastSpike[p] = t;
                                    PG.append([p,t]);
                                     
                                    cond_pre = preSynConn == p;
                                    cellsConnectedFromP = np.extract(cond_pre, postSynConn);
                                    delaysConnectedFromP = np.extract(cond_pre, delayList);
                                    weightsConnectedFromP = np.extract(cond_pre, weightList);
                                    for cell_i in range(len(cellsConnectedFromP)):
                                        timing = t+delaysConnectedFromP[cell_i];
                                        if timing>=t_max:
                                            continue;
                                        post_cell = cellsConnectedFromP[cell_i];
                                        dPSP = weightsConnectedFromP[cell_i]
                                        PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+dPSP;
                                
                                
                #print PG;                   
                if  len(PG)>sPicked:
                    if plotFigure:
                        plt.clf();
                        if drawArrow:
                            plt.hold(True);
                            ax = plt.axes();
                                    
                            for t in range(t_min,t_max):
                                for p in range(nCells*2):
                                    if (PSPRecord[p][t]==1):
                                        cond_pre = preSynConn == p;
                                        cellsConnectedFromP = np.extract(cond_pre, postSynConn);
                                        delaysConnectedFromP = np.extract(cond_pre, delayList);
                                        for cell_i in range(len(cellsConnectedFromP)):
                                            timing = t+delaysConnectedFromP[cell_i];
                                            post_cell = cellsConnectedFromP[cell_i];
                                            
                                            if (sum(PSPRecord[post_cell,timing:timing+jitter+1])>=1):
                                                dy = (post_cell-p)#*0.99;
                                                dx = (delaysConnectedFromP[cell_i])#*0.99
                                                ax.arrow(t, p,dx , dy, fc='k', ec='k', head_width=0.05, head_length=0.1)
                    
                    
                        print s_list
                        print PG;
                        PGTrans = zip(*PG);
                        plt.plot(PGTrans[1], PGTrans[0], 'o', markerfacecolor = 'w');
                        plt.ylabel('Cell Index');
                        plt.xlabel('time [ms]');
                        plt.title('supported PG '+str(s_list));
                        if PGType==1:
                            txt = "Jitter=" + str(jitter) + "ReflPeriod=" + str(reflPeriod);# + "\n" +str(PG);
                            fig.text(.1,0.01,txt)
                            svname=os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGS_Jitter="+str(jitter)+"_T="+str(s_list)+"_i="+str(i)+"len="+str(len(PG));
                        elif PGType==2:
                            txt = "spikingTh=" + str(spikingTh) + "decayRatio=" + str(decayRatio);# + "\n" +str(PG);
                            fig.text(.1,0.01,txt)
                            svname=os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGA_spikingTh="+str(spikingTh)+"_T="+str(s_list)+"_i="+str(i)+"len="+str(len(PG));                        
                        fig.savefig(svname+".png");
                    
                    PGCount+=1;
                    PGLen.append(len(PG));
                    
                    if PGType==1:
                        csvFileName = os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGS_Jitter="+str(jitter)+".txt";
                    elif PGType==2:
                        csvFileName = os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGA_spikingTh="+str(spikingTh)+".txt";
                    fd = open(csvFileName,'a')
                    fd.write(str(s_list)+" len=" +str(len(PG))+" PG="+str(PG)+"\n");
                    fd.close();
                    #plt.show();
                    
                    #plt.ylim([0, layerGDim * layerGDim - 1])
                    #plt.xlim([timeBegin, timeBegin+simulationTime])
                    
                    #tmp = input();
                        
    print PGCount;                    
    print mean(PGLen);
