from brian2 import ms,mV;
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

def runCalcPG(nStims,nTrans,analysisLayer):

    PGTypes = [3]; #1:supported PG, 2:adapted PGs, 3:activated PGs
    isTrained = False;
    plotFigure = False;    
    drawArrow = True;
    
    nCells = 100;
    sPicked = 3; 
    NbSpikesNeeded = 3;
    NbSpikesMax = 200;
    MaxTimeSpan = 100;
    reflPeriod = refractoryPeriod/ms;
    minPGLen = NbSpikesNeeded+1;

    #for alg1: PG supported
    jitter = 1;


    #for alg2: PG adopted
    triggerTh = sPicked/3;#*(gmax*0.2);#this number is arbitrary
    spikingTh = Vth_ex/mV;
    alg2_Vrev_ex = Vrev_ex/mV;
    #decayRatio = 0.5;
    alg2_V0_ex = V0_ex/mV;
    alg2_taum_ex = taum_ex/ms;
    alg2_tau_ex = tau_ex/ms;#*0.8
#     decalyRatio = (1-timeStep/taum_ex);
    w_const_FF = 1.5/gmax;#this number is arbitrary
    w_const_Lat = 1.5/conductanceConst_E2E;
    jitter_plotArrow = np.log2(0.9)/np.log2(1-1./alg2_taum_ex)
    print jitter_plotArrow;
    
    
    
    #Load Data
    analysisLayer = analysisLayer-1;
    
    f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/"+ str(trainingEpochs+1) + "_netStates_E2E.pkl", "rb");
    netState_E2E = pickle.load(f);
    f.close();
    
    if isTrained:
        fileName_netState = os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/"+ str(trainingEpochs+1) + "_netStates_L2L.pkl";
    else:
        fileName_netState = os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/0_netStates_L2L.pkl";
    fileName_spikes = os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/" + str(trainingEpochs+1)+'_spikes_e.pkl';
        
    f = open(fileName_netState, "rb");
    netState_L2L = pickle.load(f);
    f.close();
    
    f = open(fileName_spikes, "rb");
    spikes_e = pickle.load(f);
    f.close();
    
    
    
    preSynConn_FF_list = [];
    postSynConn_FF_list = [];
    weightList_FF_list = [];
    delayList_FF_list = [];
    
    for l in range(nLayers-1):
        preSynConn_FF_list.append(netState_L2L[l]['i_pre']+nCells*l);
        postSynConn_FF_list.append(netState_L2L[l]['i_post']+nCells*(1+l));
        weightList_FF_list.append(netState_L2L[l]['w']);
        delayList_FF_list.append((netState_L2L[l]['delay']/ms).astype(int));

    preSynConn_RC_list = [];
    postSynConn_RC_list = [];
    weightList_RC_list = [];
    delayList_RC_list = [];

    for l in range(nLayers):
        preSynConn_RC_list.append(netState_E2E[l]['i_pre']+nCells*l);
        postSynConn_RC_list.append(netState_E2E[l]['i_post']+nCells*l);
        weightList_RC_list.append(netState_E2E[l]['w']);
        delayList_RC_list.append((netState_E2E[l]['delay']/ms).astype(int));
    
    #net state of the final output layer
#     preSynConn_RC_o = netState_E2E[analysisLayer-1]['i_pre']+nCells;
#     postSynConn_RC_o = netState_E2E[analysisLayer-1]['i_post']+nCells;
#     weightList_RC_o = netState_E2E[analysisLayer-1]['w'];
#     delayList_RC_o = (netState_E2E[analysisLayer-1]['delay']/ms).astype(int);
    
    preSynConn = np.concatenate((np.concatenate(preSynConn_FF_list),np.concatenate(preSynConn_RC_list)));
    postSynConn = np.concatenate((np.concatenate(postSynConn_FF_list),np.concatenate(postSynConn_RC_list)));
    weightList = np.concatenate((np.concatenate(weightList_FF_list),np.concatenate(weightList_RC_list)));
    delayList = np.concatenate((np.concatenate(delayList_FF_list),np.concatenate(delayList_RC_list)));
    
    #spikeTime =  spikes_e[analysisLayer-2];
    
    #init
    PGLen = [];
    PGTriggers = [];
    PGList = [];
    
    
    if plotFigure:
        fig = plt.figure();
        

    for PGType in PGTypes:
        if PGType==1 or PGType==2:
            if PGType==1:
                fileName = os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGS_Jitter="+str(jitter);
            elif PGType==2:
                fileName = os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGA_Trained="+ str(isTrained)+"_spikingTh="+str(spikingTh);
            
#             for s_list in [(2, 3, 50)]:#list(itertools.combinations(range(nCells),sPicked)):
            for s_list in list(itertools.combinations(range(nCells),sPicked)):
    #             if len(PGList)>10:
    #                 break;
                #print s_list
                #look for PGs triggered by this combination
                print s_list
                        
                
                triggeringConnections = [];
                for i in range(nCells*2):
                    triggeringConnections.append([]);
                if PGType == 2:
                    gSumList=np.zeros(nCells*2);
                    
                
                #count connections comming from triggers, to find common triggers output neurons
                for p in s_list:
                    cond_pre = preSynConn_FF_list[0] == p;
                    listConnectedFromTp = np.extract(cond_pre, postSynConn_FF_list[0]);
                    for i in listConnectedFromTp:
                        if p not in triggeringConnections[i]:#remove duplicate
                            triggeringConnections[i].append(p);
                            if PGType==2:
                                cond_post = postSynConn_FF_list[0] == i;
                                weightListPtoI = np.extract(cond_pre & cond_post, weightList_FF_list[0]);
                                gSumList[i]=gSumList[i]+weightListPtoI[0]*w_const_FF;
                
                
                for i in range(nCells,nCells*2): 
                    PG = []
                    if PGType==1:
                        PSPTable = np.zeros((nCells*nLayers,MaxTimeSpan));
                    else:
                        PSPTable = np.zeros((nCells*nLayers,MaxTimeSpan))+alg2_V0_ex
                        geTable = np.zeros((nCells*nLayers,MaxTimeSpan));
                    PSPRecord = np.zeros((nCells*nLayers,MaxTimeSpan));
                    tLastSpike = np.zeros(nCells*nLayers)-reflPeriod;
                    
                    #if there is satisfactory incoming connections from the trigger..
                    if ((len(triggeringConnections[i]) >= sPicked) and (PGType==1 or gSumList[i]>triggerTh)):
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
                            if PGType==1:
                                for cell_i in range(len(cellsConnectedFromP)):
#                                     if p is not cellsConnectedFromP[cell_i]:                                    
                                    timing = initSpikeTimings[p]+delaysConnectedFromP[cell_i];
                                    post_cell = cellsConnectedFromP[cell_i];
                                    PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+1;
                            elif PGType==2:
                                weightsConnectedFromP = np.extract(cond_pre, weightList);
                                for cell_i in range(len(cellsConnectedFromP)):
                                    timing = initSpikeTimings[p]+delaysConnectedFromP[cell_i];
                                    post_cell = cellsConnectedFromP[cell_i];
                                    ge = weightsConnectedFromP[cell_i]
                                    
                                    layer_pre = int(np.floor((p*1.0)/nCells));
                                    layer_post = int(np.floor((cellsConnectedFromP[cell_i]*1.0)/nCells));
                                                                         
                                    if layer_pre == layer_post: 
                                        geTable[post_cell,timing]= geTable[post_cell,timing]+ge*w_const_Lat;
                                    else:
                                        geTable[post_cell,timing]= geTable[post_cell,timing]+ge*w_const_FF;
                                        #PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+dv;
                                        #print(PSPTable[post_cell,timing]);
                                        #PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+dPSP;
                         
                        for t in range(MaxTimeSpan):
                            if len(PG)<NbSpikesMax:
                                for p in range(nCells*nLayers):
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
                                                if timing>=MaxTimeSpan:# or p == cellsConnectedFromP[cell_i]:
                                                    continue;
                                                post_cell = cellsConnectedFromP[cell_i];
                                                PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+1;
                                    elif PGType==2:
                                        if t>0:
                                            #PSPTable[p,t]=PSPTable[p,t-1]*decayRatio+PSPTable[p,t];
                                            #dge/dt = -ge/tau_ex
                                            geTable[p,t] = geTable[p,t] + geTable[p,t-1]*(1-(1./alg2_tau_ex))
                                            dv = ((alg2_V0_ex -PSPTable[p,t-1]) + (geTable[p,t] * (alg2_Vrev_ex-PSPTable[p,t-1])))/ alg2_taum_ex ;
                                            PSPTable[p,t]=PSPTable[p,t-1] + dv; 
                                            
                                            
                                            if ((PSPTable[p,t]>=spikingTh) and (t>=tLastSpike[p]+reflPeriod)):
                                                PSPRecord[p,t]=1;
                                                PSPTable[p,t] = alg2_V0_ex;#reset PSP
                                                tLastSpike[p] = t;
                                                PG.append([p,t]);
                                                 
                                                cond_pre = preSynConn == p;
                                                cellsConnectedFromP = np.extract(cond_pre, postSynConn);
                                                delaysConnectedFromP = np.extract(cond_pre, delayList);
                                                weightsConnectedFromP = np.extract(cond_pre, weightList);
                                                for cell_i in range(len(cellsConnectedFromP)):
                                                    timing = t+delaysConnectedFromP[cell_i];
                                                    if timing>=MaxTimeSpan or p == cellsConnectedFromP[cell_i]:
                                                        continue;
                                                    post_cell = cellsConnectedFromP[cell_i];
                                                    ge = weightsConnectedFromP[cell_i]
                                                    #PSPTable[post_cell,timing]=PSPTable[post_cell,timing]+ge;
                                                    
                                                    layer_pre = int(np.floor((p*1.0)/nCells));
                                                    layer_post = int(np.floor((post_cell*1.0)/nCells));
                                                    
                                                    if layer_post==layer_pre:
                                                        geTable[post_cell,timing]=geTable[post_cell,timing]+ge*w_const_Lat;
                                                    else:
                                                        geTable[post_cell,timing]=geTable[post_cell,timing]+ge*w_const_FF;

                                                    #PSPTable[post_cell,timing]=((alg2_Vrev_ex-PSPTable[post_cell,timing])*ge)/alg2_taum_ex;
                                            
                        #print nad save PG;                   
                        if  len(PG)>=minPGLen:#sPicked:
                            if plotFigure and len(PG)>20 and len(PG)<150:
                                plt.clf();
                                if drawArrow:
                                    plt.hold(True);
                                    ax = plt.axes();
                                            
                                    for t in range(MaxTimeSpan):
                                        for p in range(nCells*nLayers):
                                            if (PSPRecord[p][t]==1):
                                                cond_pre = preSynConn == p;
                                                cellsConnectedFromP = np.extract(cond_pre, postSynConn);
                                                delaysConnectedFromP = np.extract(cond_pre, delayList);
                                                for cell_i in range(len(cellsConnectedFromP)):
                                                    timing = t+delaysConnectedFromP[cell_i];
                                                    post_cell = cellsConnectedFromP[cell_i];
                                                    if PGType==1:
                                                        if (sum(PSPRecord[post_cell,timing:timing+jitter+1])>=1):
                                                            dy = (post_cell-p)#*0.99;
                                                            dx = (delaysConnectedFromP[cell_i])#*0.99
                                                            ax.arrow(t, p,dx , dy, fc='k', ec='k', head_width=0.05, head_length=0.1)
                                                    elif PGType==2:
                                                        if (sum(PSPRecord[post_cell,timing:timing+int(jitter_plotArrow)+1])>=1):
                                                            dy = (post_cell-p)#*0.99;
                                                            dx = (delaysConnectedFromP[cell_i])#*0.99
                                                            ax.arrow(t, p,dx , dy, fc='k', ec='k', head_width=0.05, head_length=0.1)
                                    
                                    for l in range(1,nLayers):
                                        plt.plot((0,np.max(tLastSpike)*1.1),(nCells*l,nCells*l),'k--');
                                    plt.xlim((0,np.max(tLastSpike)*1.1))
                            
                            
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
                                    txt = "spikingTh=" + str(spikingTh) + "alg2_taum_ex=" + str(alg2_taum_ex);# + "\n" +str(PG);
                                    fig.text(.1,0.01,txt)
                                    svname=os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGA_Trained="+ str(isTrained) +"_spikingTh="+str(spikingTh)+"_T="+str(s_list)+"_i="+str(i)+"len="+str(len(PG));                        
                                fig.savefig(svname+".png");
                            
                            PGLen.append(len(PG));
                            PGTriggers.append([len(PG),PG[0],PG[1],PG[2],PG[3]]);
                            PGList.append(PG);
                            
                            
                            fd = open(fileName+".txt",'a')
                            fd.write(str(s_list)+" len=" +str(len(PG))+" PG="+str(PG)+"\n");
                            fd.close();
                            #plt.show();
                            
                            #plt.ylim([0, layerGDim * layerGDim - 1])
                            #plt.xlim([timeBegin, timeBegin+simulationTime])
                            
                            #tmp = input();
              
                          
            print len(PGList);                    
            print np.mean(PGLen);
                
            fd = open(fileName+"_Triggers.pkl", "wb");
            pickle.dump(PGTriggers, fd);
            fd.close();
            
            fd = open(fileName+"_PGList.pkl", "wb");
            pickle.dump(PGList, fd);
            fd.close();
            
            fd = open(fileName+".txt",'a')
            fd.write("\n PGCount=" +str(len(PGList))+" AvgPGLen="+str(np.mean(PGLen))+"\n");
            fd.close();
            
            
            
            
            
        if PGType == 3:
    #         #load triggers
    #         fileName = os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGA_Trained="+ str(isTrained)+"_spikingTh="+str(spikingTh)+"_Triggers.pkl";
    #         f = open(fileName, "rb");
    #         PGTriggers = pickle.load(f);
    #         f.close();
    #
            fileName = os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGA_Trained="+ str(isTrained)+"_spikingTh="+str(spikingTh)+"_PGList.pkl";
            f = open(fileName, "rb");
            PGList = pickle.load(f);
            f.close();
            
            triggerActivatedCount = np.zeros((nStims,len(PGList)));
            
            #init
            phases = []
            if isTrained:
                for stim in range(nStims):
                    phases.append((testingTime+trainingTime*trainingEpochs)*nStims*nTrans+stim*testingTime*nTrans);
            else:
                for stim in range(nStims):
                    phases.append(testingTime*stim*nTrans);
            
    
            phase_i = 0;
    #         triggerList = [];
            for t_min in phases:
                print str(t_min)
                phase_i = phase_i+1;
                triggeredList = []
                PGLen_sum = 0;
                t_max = t_min+testingTime*nTrans;
                PG_index = 0;
                print len(PGList);
                for PG in PGList:
                    triggeringPGLen = len(PG);#len of triggering PG
                    triggeredCell = PG[sPicked];#first cell triggered by the triggers
                    triggers=np.array(PG[0:sPicked]);
                    triggers_sorted= triggers[np.argsort(triggers[:,1])];
                    p_init = triggers_sorted[0][0];
                    t_init = triggers_sorted[0][1];
                        
                    spikeTime_pre =  spikes_e[0][p_init];
                    cond_pre = (t_min < spikeTime_pre) & (spikeTime_pre < t_max);
                    spikeTime_firstTrigger=np.extract(cond_pre, spikeTime_pre);    #extract spike timings between specified timing
    
                    for t in spikeTime_firstTrigger:#for all spikes emmitted by the first trigger
                        areTriggered = True;
                        for Np in triggers_sorted[1:sPicked]:
                            cell_i = int(Np[0]%nCells);
                            layer = int(np.floor(Np[0]*1.0/nCells));
                            
#                             print "pre: " + str(layer) + " " +str(cell_i);
                            
                            spikeTime_p = spikes_e[layer][cell_i];
                            cond =  (t+Np[1]-jitter <= spikeTime_p) & (spikeTime_p <= t+Np[1]+jitter)
                            if len(np.extract(cond, spikeTime_p))<1:
                                areTriggered = False;
                                break;
                        
                        
                        cell_i = int(triggeredCell[0]%nCells);
                        layer = int(np.floor(triggeredCell[0]*1.0/nCells))
#                         print "post: " + str(layer) + " " +str(cell_i);
                        spikeTime_post =  spikes_e[layer][cell_i];
                        cond =  (t+triggeredCell[1]-jitter <= spikeTime_post) & (spikeTime_post <= t+triggeredCell[1]+jitter);
                        if len(np.extract(cond, spikeTime_post))<1:
                            areTriggered = False;
                        
                        if areTriggered:
                            PGLen_sum=PGLen_sum+triggeringPGLen;
                            triggerActivatedCount[phase_i-1,PG_index]=triggerActivatedCount[phase_i-1,PG_index]+1;
                            triggeredList.append([]);
                            for Np in triggers_sorted:
                                triggeredList[len(triggeredList)-1].append([Np[0],t+Np[1]]);
    #             triggerList.append(PG);
                    PG_index = PG_index+1;


                corrcoef = np.corrcoef(triggerActivatedCount);
                
                
                print str(phase_i) + ": PGLen=" + str(PGLen_sum)# + " PGList=" + str(PG);
                fileName = os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PGAc_Trained="+ str(isTrained)+"_spikingTh="+str(spikingTh);
                fd = open(fileName+".txt",'a')
                fd.write(str(phase_i) + ": PGLen_sum=" + str(PGLen_sum)+"\n");# + " TriggeredList=" + str(triggeredList) + "\n");
#                 fd.write(str(phase_i) + ": PGLen_sum=" + str(PGLen_sum)+" TriggeredList=" + str(triggeredList) + "\n");
                fd.close(); 
                
                fd = open(fileName+"_PGList.pkl",'wb');
                pickle.dump(triggeredList, fd);
                fd.close();
            
#             print "triggerActivatedCount=" + str(triggerActivatedCount);
#             print "corrcoef=" + str(corrcoef); 
#             fd = open(fileName+".txt",'a')
#             fd.write("triggerActivatedCount=" + str(triggerActivatedCount)+"\n corrcoef=" + str(corrcoef));# + " TriggeredList=" + str(triggeredList) + "\n");
#             fd.close();
#             fig = plt.figure();
#             plt.imshow(corrcoef, cmap='gray', vmin=0, vmax=1, interpolation='none')
#             fig.savefig(fileName+".png");
                
               


