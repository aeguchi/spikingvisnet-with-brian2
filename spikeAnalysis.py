from brian2 import ms;
import numpy as np;
import pylab as plt;
import pickle
import os;
import errno
import InfoAnalysis

#infoAnalysis?
#ia = InfoAnalysis.InfoAnalysis(globals())
#ia.singleCellInfoAnalysis(phases = ['FR_0_blank.pkl']);
#ia.singleCellInfoAnalysis(numBins=3,weightedAnalysis=1);
#plt.show();


#nBins = 10;
# analysing cell index==0 for test development

#t_min = 0#840000#15000;#319000;
#t_max = 20000#860000#320000;

def loadParams(borrowed_globals):
    globals().update(borrowed_globals);




def runSpikeAnalysis(nStims,nTrans,analysisLayer, PIcalcOn=True,polyAnalysisOn = False,polyHist = False):
    def traceDelay(poly_indexs,poly_delays,index,delay):
        #print str(index) + " " + str(delay);
        polyTableConnected = polyTable[index][preList];
        preMax = polyTableConnected.max();
        #print "max: " + str(preMax);
        if preMax<polyChainDetectTh:
            return;
        preArgMax = np.where(polyTableConnected==preMax);
        
        for in_max in range(len(preArgMax[0])):
            cond1 = np.array(poly_indexs)==preList[preArgMax[0][in_max]];#check if the pair of the values already exists
            cond2 = np.array(poly_delays)==delay+int(preArgMax[1][in_max]);
            cond = cond1 & cond2;
            #print str(preArgMax[0][in_max]) + ", " +str(delay+preArgMax[1][in_max]);
            #if preArgMax[1][in_max]!=0 and delay+preArgMax[1][in_max]<maxDelay:
            if not (True in cond) and delay+preArgMax[1][in_max]<maxDelay:
                print str(preList[preArgMax[0][in_max]]) + " " + str(delay+int(preArgMax[1][in_max]));
                poly_indexs.append(preList[preArgMax[0][in_max]]);
                poly_delays.append(delay+int(preArgMax[1][in_max]));
                traceDelay(poly_indexs,poly_delays,preList[preArgMax[0][in_max]],delay+int(preArgMax[1][in_max]));
    
    
    
    
    #params
    polyChainDetectTh = 40;
    nTopCellsChosenForPI = 10#10;
    nCells = layerDim*layerDim;
    maxDelay = delayConst_connBottomUp/ms;
    nSpikesUsedForPI = maxDelay*testingTime*nTrans/1000
    #experimentName = 'BO_single1_gabMod2';
    
    
    #init
    SpikeData = [str(trainingEpochs+1)+'_spikes_e.pkl', str(trainingEpochs+1)+'_spikes_i.pkl', str(trainingEpochs+1)+'_pikes_b.pkl', str(trainingEpochs+1)+'_pikes_g.pkl'];
    spikes_e = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/" + SpikeData[0], "rb"));
    netState_L2L = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/"+ str(trainingEpochs+1) + "_netStates_L2L.pkl", "rb"));
    
    
    #analysisLayer = nLayers-1;#(3)
    analysisLayer = analysisLayer-1;
    
    preSynConn = netState_L2L[analysisLayer-1]['i_pre'];
    postSynConn = netState_L2L[analysisLayer-1]['i_post'];
    weightList = netState_L2L[analysisLayer-1]['w'];
    delayList = netState_L2L[analysisLayer-1]['delay'];
    
    phases = []
    PI = np.zeros((2,nStims,nCells));#phase(blank, trained), number of stimuli, number of cells
    
    for stim in range(nStims):
        phases.append(testingTime*stim*nTrans);
    
    for stim in range(nStims):
        phases.append((testingTime+trainingTime*trainingEpochs)*nStims*nTrans+stim*testingTime*nTrans);
    
    PI_phase = 0;
    stim_index = 0;
    for t_min in phases:
        stim_index = stim_index%nStims;#param needed to save PI
        t_max = t_min+testingTime*nTrans;
        if t_min>testingTime*nStims*nTrans:
            PI_phase = 1;#param needed to save PI
        else:
            PI_phase = 0;
        
        #fig_hist = plt.figure(1 , figsize=(40, 40),dpi=500);
        #fig_poly = plt.figure(2 , figsize=(40, 40),dpi=100);
        
        #calculate FR in layer 0:
        meanFR_pre = 0;
        meanFR_post = 0;
        for i in range(100):
            spikeTime_pre = spikes_e[analysisLayer-1][i]; #spiking time of cell i_post in post synaptic layer (layer 1)
            cond_post1 = spikeTime_pre < t_max;
            cond_post2 = t_min < spikeTime_pre;
            cond_post = cond_post1 & cond_post2;
            spikeTime_pre_ext=np.extract(cond_post, spikeTime_pre)    #extract spike timings between specified timing
            meanFR_pre += len(spikeTime_pre_ext)*1.0/nCells/(t_max-t_min)*1000;
            
            
            spikeTime_post = spikes_e[analysisLayer][i]; #spiking time of cell i_post in post synaptic layer (layer 1)
            cond_post1 = spikeTime_post < t_max;
            cond_post2 = t_min < spikeTime_post;
            cond_post = cond_post1 & cond_post2;
            spikeTime_post_ext=np.extract(cond_post, spikeTime_post)    #extract spike timings between specified timing
            meanFR_post += len(spikeTime_post_ext)*1.0/nCells/(t_max-t_min)*1000;
        print meanFR_pre;
        print meanFR_post;
        
        
        for i_post in range(100):
            print i_post
            #plt.clf();
            #prevSpikeTime = np.ones([len(preList), len(spikes_e[1][i_post])])*(maxDelay+1);  # pre cell index, number of spikes
            prevSpikeTime = np.ones([nCells, len(spikes_e[analysisLayer][i_post])])*(maxDelay+1);  # pre cell index, number of spikes
            relativeSpikeTime = []; #nCells,
            
            spikeTime_post = spikes_e[analysisLayer][i_post]; #spiking time of cell i_post in post synaptic layer (layer 1)
            cond_post1 = spikeTime_post < t_max;
            cond_post2 = t_min < spikeTime_post;
            cond_post = cond_post1 & cond_post2;
            spikeTime_post_ext=np.extract(cond_post, spikeTime_post)    #extract spike timings between specified timing    
            
            if(len(spikeTime_post_ext)>0):  #if post synaptic cell ever spikes
                for i_pre in range(100):#for each presynaptic cell
                    count = 0; #keep track with spike index
                    relativeSpikeTime.append([]);#append cell index    
                    for t_post in spikeTime_post_ext: #for each spikes of the post synaptic cell (i_post)
                        relativeSpikeTime[i_pre].append([]);
                        if len(spikes_e[analysisLayer-1][i_pre]) > 0: #if the presynaptic cell in layer 0 ever spikes
                            # print i_pre
                            spikeTime_pre = spikes_e[analysisLayer-1][i_pre]; #store the spike timings of the presynaptic cell in layer 0
                            spikeTime_pre = spikeTime_pre[::-1];    #and reverse order
                            cond_pre1 = t_post-maxDelay < spikeTime_pre;
                            cond_pre2 = spikeTime_pre < t_post;
                            cond = cond_pre1 & cond_pre2;   #extract presynaptic spike timing between the timing of post spike and XX ms before the spike
                            spikeTime_pre_ext = np.extract(cond, spikeTime_pre);
                            if (len(spikeTime_pre_ext) > 0):
                                t_back = t_post - spikeTime_pre_ext[0];#t_post - last activation of the presynaptic cells within the interval
                                prevSpikeTime[i_pre][count] = t_back;
                                for spike_i in range(len(spikeTime_pre_ext)):
                                    relativeSpikeTime[i_pre][count].append(t_post - spikeTime_pre_ext[spike_i]);
                        count += 1;
                    #print count
                
         
        
                #then, try to trace
                if polyHist:
                    
                    try:
                        os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/hist"+str(t_min)+"-"+str(t_max));
                    except OSError as exception:
                        if exception.errno != errno.EEXIST:
                            raise
        
                    
                    fig_hist = plt.figure(1 , figsize=(40, 40),dpi=500);
                    plt.clf();
                    for i_pre in range(100):
                        plt.subplot(10,10,i_pre+1);
                        c = 'b';
                        #overlay actual connectivity
                        cond_ipost = postSynConn == i_post
                        cond_ipre = preSynConn == i_pre
                        cond_syn = cond_ipost&cond_ipre
                        preList = np.extract(cond_syn,preSynConn);
                        weight = np.extract(cond_syn,weightList);
                        delay = np.extract(cond_syn, delayList)/ms;
                        if len(preList)>0:
                            plt.title('connected w:' + "{:.2f}".format(weight[0]/gmax));
                            c = 'r';
            
                        tmp = np.extract(prevSpikeTime[i_pre]<maxDelay+1,prevSpikeTime[i_pre]);
            
                        if len(tmp)>1:
                            plt.hist(tmp,bins=10,range = (0,maxDelay),color=c);
                        elif len(tmp)==1:
                            plt.bar(tmp,1,align="center");
                            
                        if len(preList)>0:
                            plt.hold(True);
                            plt.plot(delay[0],len(spikeTime_post_ext)/2,'*');
                            
                        plt.ylim(len(spikeTime_post_ext)*-0.1,len(spikeTime_post_ext)*1.1)
                        plt.xlim(-0.01*maxDelay,maxDelay*1.01);
                        
                        #plt.show()
                        
                    fig_hist.savefig(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/hist"+str(t_min)+"-"+str(t_max)+"/hist_"+str(t_min)+"-"+str(t_max)+"_"+str(i_post)+".png");
        
        
                if PIcalcOn:
                    #calculating PI
                    if count>=nSpikesUsedForPI:
                        polyTableTmp = np.zeros([nCells,maxDelay+1]);
                        randIndex = np.random.permutation(count);
                        cond_ipost = postSynConn == i_post
                        preList = np.extract(cond_ipost,preSynConn);
                        
                        for i_SpikeTrain in randIndex[0:nSpikesUsedForPI]:
                            #print len(randIndex[0:nSpikesUsedForPI]);
                            for inputCell in range(100):
                                for spikeTime in relativeSpikeTime[inputCell][i_SpikeTrain]:
                                   polyTableTmp[inputCell][int(spikeTime)] += 1; 

                        sortedPolyTable = np.sort(polyTableTmp.reshape((maxDelay+1)*nCells));
                        if meanFR_pre>0:
                            PI[PI_phase,stim_index,i_post]=np.mean(sortedPolyTable[-1-nTopCellsChosenForPI:-1])/(meanFR_pre*nSpikesUsedForPI/1000);
        
                
                #for each input spikes find the correlation
                #polyTable = np.zeros([nCells,nCells,maxDelay]);#post,pre,diff
                if polyAnalysisOn:
                    try:
                        os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/poly"+str(t_min)+"-"+str(t_max));
                    except OSError as exception:
                        if exception.errno != errno.EEXIST:
                            raise
                        
                    fig_poly = plt.figure(2 , figsize=(20, 20),dpi=100);
                    plt.clf();
                    plt.title("Post Synaptic Cell index: "+str(i_post));
                    
                    cond_ipost = postSynConn == i_post
                    preList = np.extract(cond_ipost,preSynConn);
                    preList_sorted = np.sort(preList);
                    polyTable = np.zeros([nCells,nCells,maxDelay+1]);#post,pre,diff
                    spikeCount = np.zeros(nCells);
                    polyTableTmp = np.zeros([nCells,maxDelay+1]);
                    subplotDim = np.ceil(np.sqrt(len(preList)+1));
                    
                    
                    if count>0:
                        for i_SpikeTrain in range(count):#for each post synaptic spike
                            #print relativeSpikeTime[0][i_SpikeTrain]
                            for inputCell in range(100):
                                for spikeTime in relativeSpikeTime[inputCell][i_SpikeTrain]:
                                   polyTableTmp[inputCell][int(spikeTime)] += 1; 

                        #plt.subplot(int(subplotDim)+1,int(subplotDim),1);
                        polyTableTmpConnected = polyTableTmp[preList_sorted];
                        plt.subplot(int(subplotDim)+1,2,1);
                        plt.imshow(polyTableTmpConnected,interpolation='none');                
                        plt.colorbar()
                        plt.title("postSynapticCell:" + str(i_post) + " count:" + str(count));
                        #plt.title("(PostSynCell: " + str(i_post) + ") Distribution of the delay of spikes after spikes of PreSynCell " + str(i))
                        plt.xlabel("delay [ms]");
                        #plt.ylabel("index of neuron in layer 0");
                        plt.yticks(range(len(preList)), preList_sorted)
                        plt.gca().invert_xaxis()
                        
                        
                        subplot_i=0;
                        for cell_focus in preList:
                            subplot_i+=1;
                            spikeTime_focus = spikes_e[analysisLayer-1][cell_focus]; #spiking time of cell i_post in post synaptic layer (layer 1)
                            cond_focus1 = spikeTime_focus < t_max;
                            cond_focus2 = t_min < spikeTime_focus;
                            cond_focus = cond_focus1 & cond_focus2;
                            spikeTime_focus_ext=np.extract(cond_focus, spikeTime_focus)
                            for i_pre in preList:    
                                for t_focus in spikeTime_focus_ext: #for each spikes of the post synaptic cell (i_post)
                                    if len(spikes_e[analysisLayer-1][i_pre]) > 0: #if the presynaptic cell in layer 0 ever spikes
                                        # print i_pre
                                        spikeTime_pre = spikes_e[analysisLayer-1][i_pre]; #store the spike timings of the presynaptic cell in layer 0
                                        spikeTime_pre = spikeTime_pre[::-1];    #and reverse order
                                        cond_pre1 = t_focus-maxDelay < spikeTime_pre;
                                        cond_pre2 = spikeTime_pre < t_focus;
                                        cond = cond_pre1 & cond_pre2;   #extract presynaptic spike timing between the timing of post spike and XX ms before the spike
                                        spikeTime_pre_ext = np.extract(cond, spikeTime_pre);
                                        if (len(spikeTime_pre_ext) > 0):
                                            for spike_i in range(len(spikeTime_pre_ext)):
                                                t_back = int(t_focus - spikeTime_pre_ext[spike_i]);
                                                polyTable[cell_focus][i_pre][t_back] += 1;
              
                            plt.subplot(int(subplotDim)+1,int(subplotDim),subplot_i+int(subplotDim));

                            plt.imshow(polyTable[cell_focus][preList_sorted],interpolation='none');
                            plt.colorbar()
                            plt.title("i:" + str(cell_focus) + " max:" + str(spikeCount[cell_focus]));
                            plt.xlabel("delay [ms]");
                            plt.gca().invert_xaxis()
                            plt.yticks(range(len(preList)), preList_sorted)
                        
                        
                        #tracing poly
                        poly_indexs = []
                        poly_delays = []
                        
                        postMax = polyTableTmpConnected.max();
                        postArgMax = np.where(polyTableTmpConnected==postMax);
                        
                        
                        if postMax>=polyChainDetectTh:
                            for max_in in range(len(postArgMax[0])):
                                poly_indexs.append(preList_sorted[postArgMax[0][max_in]]);
                                poly_delays.append(int(postArgMax[1][max_in]));
                                traceDelay(poly_indexs,poly_delays,preList_sorted[postArgMax[0][max_in]],int(postArgMax[1][max_in]))
                                plt.subplot(int(subplotDim)+1,2,2);
                                plt.plot(poly_delays,poly_indexs,'*');
                                plt.xlim([maxDelay*-0.01,maxDelay*1.01]);
                                plt.ylim([-1,nCells]);
                                plt.gca().invert_yaxis()
                                plt.gca().invert_xaxis()
                                #plt.show()
                         
                        fig_poly.savefig(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/poly"+str(t_min)+"-"+str(t_max)+"/poly_"+str(t_min)+"-"+str(t_max)+"_"+str(i_post)+"_polyLen"+str(len(poly_indexs))+".png");#"_PI{:.2f}".format(PI[i_post])+".png");
            
                        
            
            #print "stop";
            #plt.plot(prevSpikeTime[:][1::]);
            
            #plt.show()
        
        if PIcalcOn:
            PImean = np.mean(np.extract(PI[PI_phase][stim_index]!=0,PI[PI_phase][stim_index]));
            print PImean;
            f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PI_"+str(t_min)+"-"+str(t_max)+".txt","w");    
            f.write(str(PImean));
            f.close();
        
        # #spikes[layer][cellIndex]
        # layer = 1;
        # spikes = tmp[layer];
        # endingTime = 7200; 
        # print endingTime
        # plt.subplot(3,1,1);
        # plt.hold(True)
        # for i in range(100):
        #     plt.plot(spikes[i]/ms,np.ones(len(spikes[i]))*i,'.');
        # plt.xlim([endingTime-3000,endingTime-2000])
        # 
        # plt.subplot(3,1,2);
        # plt.hold(True)
        # for i in range(100):
        #     plt.plot(spikes[i]/ms,np.ones(len(spikes[i]))*i,'.');
        # plt.xlim([endingTime-2000,endingTime-1000])
        # 
        # plt.subplot(3,1,3);
        # plt.hold(True)
        # for i in range(100):
        #     plt.plot(spikes[i]/ms,np.ones(len(spikes[i]))*i,'.');
        # plt.xlim([endingTime-1000,endingTime])
        # 
        # plt.show()
        stim_index+=1;
        
    
    if PIcalcOn:
        PImeanBlank = 0;
        PImeanTrained = 0;
        for stim in range(nStims):
            tmp = np.mean(np.extract(PI[0][stim]!=0,PI[0][stim]));
            if np.isnan(tmp) == False:
                PImeanBlank += tmp;
            tmp = np.mean(np.extract(PI[1][stim]!=0,PI[1][stim]));
            if np.isnan(tmp) == False:
                PImeanTrained += tmp;
        
        print "PImeanImprovement: " + str(PImeanTrained-PImeanBlank) + " (blank:" + str(PImeanBlank/nStims) + ", trained:" + str(PImeanTrained/nStims) + ")";
        f = open(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/PI_improvement.txt","w");    
        f.write(str(PImeanBlank-PImeanTrained));
        f.close();
        
        f = open(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/infoAndPI.txt","a");
        f.write(str(PImeanBlank-PImeanTrained)+" ");
        f.close();
