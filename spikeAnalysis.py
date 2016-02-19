from brian2 import ms;
import numpy as np;
import pylab as plt;
import pickle
import os;
import errno
import InfoAnalysis

experimentName = 'BO_single1_1000';

#infoAnalysis?
ia = InfoAnalysis.InfoAnalysis(globals())
ia.singleCellInfoAnalysis(phases = ['FR_0_blank.pkl']);
#ia.singleCellInfoAnalysis();
plt.show();




#SpikeData = ['0_spikes_e.pkl', '0_spikes_i.pkl', '0_pikes_b.pkl', '0_pikes_g.pkl'];
SpikeData = ['101_spikes_e.pkl', '101_spikes_i.pkl', '101_pikes_b.pkl', '101_pikes_g.pkl'];
spikes_e = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/" + SpikeData[0], "rb"));
#netState_L2L = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/0_netStates_L2L.pkl", "rb"));
netState_L2L = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/101_netStates_L2L.pkl", "rb"));

preSynConn = netState_L2L[0]['i_pre'];
postSynConn = netState_L2L[0]['i_post'];
weightList = netState_L2L[0]['w'];
delayList = netState_L2L[0]['delay'];

polyAnalysisOn = True
polyHist = False

#nBins = 10;
# analysing cell index==0 for test development
nCells = 100;
t_min = 90000#15000;#319000;
t_max = 100000#320000;
max_delay = 20;
gmax = 2#2.5
polyChainDetectTh = 10;


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
        #if preArgMax[1][in_max]!=0 and delay+preArgMax[1][in_max]<max_delay:
        if not (True in cond) and delay+preArgMax[1][in_max]<max_delay:
            if preList[preArgMax[0][in_max]]==58 and delay+int(preArgMax[1][in_max])==19:
                print "test"
            print str(preList[preArgMax[0][in_max]]) + " " + str(delay+int(preArgMax[1][in_max]));
            poly_indexs.append(preList[preArgMax[0][in_max]]);
            poly_delays.append(delay+int(preArgMax[1][in_max]));
            traceDelay(poly_indexs,poly_delays,preList[preArgMax[0][in_max]],delay+int(preArgMax[1][in_max]));









#fig_hist = plt.figure(1 , figsize=(40, 40),dpi=500);
#fig_poly = plt.figure(2 , figsize=(40, 40),dpi=100);

for i_post in range(100):
    print i_post
    #plt.clf();
    #prevSpikeTime = np.ones([len(preList), len(spikes_e[1][i_post])])*(max_delay+1);  # pre cell index, number of spikes
    prevSpikeTime = np.ones([nCells, len(spikes_e[1][i_post])])*(max_delay+1);  # pre cell index, number of spikes
    spikeTime_post = spikes_e[1][i_post]; #spiking time of cell i_post in post synaptic layer (layer 1)
    cond_post1 = spikeTime_post < t_max;
    cond_post2 = t_min < spikeTime_post;
    cond_post = cond_post1 & cond_post2;
    spikeTime_post_ext=np.extract(cond_post, spikeTime_post)    #extract spike timings between specified timing
    
    if(len(spikeTime_post_ext)>0):  #if post synaptic cell ever spikes
        count = 0;
        for t_post in spikeTime_post_ext: #for each spikes of the post synaptic cell (i_post)
            for i_pre in range(100):        #for each presynaptic cell
                if len(spikes_e[0][i_pre]) > 0: #if the presynaptic cell in layer 0 ever spikes
                    # print i_pre
                    spikeTime_pre = spikes_e[0][i_pre]; #store the spike timings of the presynaptic cell in layer 0
                    spikeTime_pre = spikeTime_pre[::-1];    #and reverse order
                    cond_pre1 = t_post-max_delay < spikeTime_pre;
                    cond_pre2 = spikeTime_pre < t_post;
                    cond = cond_pre1 & cond_pre2;   #extract presynaptic spike timing between the timing of post spike and XX ms before the spike
                    spikeTime_pre_ext = np.extract(cond, spikeTime_pre);
                    if (len(spikeTime_pre_ext) > 0):
                        t_back = t_post - spikeTime_pre_ext[0];
                        prevSpikeTime[i_pre][count] = t_back;
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
    
                tmp = np.extract(prevSpikeTime[i_pre]<max_delay+1,prevSpikeTime[i_pre]);
    
                if len(tmp)>1:
                    plt.hist(tmp,bins=10,range = (0,max_delay),color=c);
                elif len(tmp)==1:
                    plt.bar(tmp,1,align="center");
                    
                if len(preList)>0:
                    plt.hold(True);
                    plt.plot(delay[0],len(spikeTime_post_ext)/2,'*');
                    
                plt.ylim(len(spikeTime_post_ext)*-0.1,len(spikeTime_post_ext)*1.1)
                plt.xlim(-0.01*max_delay,max_delay*1.01);
                
                #plt.show()
                
            fig_hist.savefig(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/hist"+str(t_min)+"-"+str(t_max)+"/hist_"+str(t_min)+"-"+str(t_max)+"_"+str(i_post)+".png");

        #for each input spikes find the correlation
        #polyTable = np.zeros([nCells,nCells,max_delay]);#post,pre,diff
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
            polyTable = np.zeros([nCells,nCells,max_delay]);#post,pre,diff
            spikeCount = np.zeros(nCells);
            polyTableTmp = np.zeros([nCells,max_delay]);
            subplotDim = np.ceil(np.sqrt(len(preList)+1));
                        
            sCount = 0;
            for i_SpikeTrain in range(count):#for each post synaptic spike
                tmp = prevSpikeTime[:,i_SpikeTrain];  #store delay of each presynaptic spike (before the ith post synaptic spike)
                sCount+=1;    
                sorted = np.sort(tmp)       #sorted the delay of each input cell (before the ith post synaptic spike)
                sorted_i = np.argsort(tmp)  #sorted index
                cond_sorted1 = sorted<max_delay+1; #eliminated all the stored delay bigger than max_delay
                ext_sorted = np.extract(cond_sorted1,sorted);
                ext_sorted_i = np.extract(cond_sorted1,sorted_i);
                #print str(cell_in) + " " + str(t_backs);
                for order_in in range(len(ext_sorted)):
                    t_back = int(ext_sorted[order_in]/1);
                    polyTableTmp[ext_sorted_i[order_in]][t_back] += 1;
            

            #plt.subplot(int(subplotDim)+1,int(subplotDim),1);
            polyTableTmpConnected = polyTableTmp[preList];
            plt.subplot(int(subplotDim)+1,2,1);
            plt.imshow(polyTableTmpConnected,interpolation='none');                
            plt.colorbar()
            plt.title("postSynapticCell:" + str(i_post) + " count:" + str(count));
            #plt.title("(PostSynCell: " + str(i_post) + ") Distribution of the delay of spikes after spikes of PreSynCell " + str(i))
            plt.xlabel("delay [ms]");
            #plt.ylabel("index of neuron in layer 0");
            plt.yticks(range(len(preList)), preList)
            plt.gca().invert_xaxis()
            
            
            subplot_i=0;
            for cell_focus in preList:
                subplot_i+=1;
                for i_SpikeTrain in range(count):#for each post synaptic spike
                    tmp = prevSpikeTime[:,i_SpikeTrain];  #store delay of each presynaptic spike (before the ith post synaptic spike)
                    if tmp[cell_focus]<max_delay+1: #
                        spikeCount[cell_focus]+=1;    
                        sorted = np.sort(tmp)       #sorted the delay of each input cell (before the ith post synaptic spike)
                        sorted_i = np.argsort(tmp)  #sorted index
                        cond_sorted1 = sorted<max_delay+1; #eliminated all the stored delay bigger than max_delay
                        cond_sorted2 = tmp[cell_focus]<sorted
                        cond_sorted = cond_sorted1&cond_sorted2;
                        ext_sorted = np.extract(cond_sorted,sorted);
                        ext_sorted_i = np.extract(cond_sorted,sorted_i);
                        t_backs=ext_sorted-tmp[cell_focus];
                        #print str(cell_in) + " " + str(t_backs);
                        for order_in in range(len(ext_sorted)):
                            t_back = int(t_backs[order_in]/1);
                            polyTable[cell_focus][ext_sorted_i[order_in]][t_back] += 1;                
                plt.subplot(int(subplotDim)+1,int(subplotDim),subplot_i+int(subplotDim));
                #print polyTable[i]
                #plt.clf()
                #plt.imshow(polyTable[i],interpolation='none');
                #plt.imshow(polyTable[cell_focus][preList],interpolation='none',vmin=0,vmax=spikeCount[cell_focus]);
                plt.imshow(polyTable[cell_focus][preList],interpolation='none');
                plt.colorbar()
                plt.title("i:" + str(cell_focus) + " max:" + str(spikeCount[cell_focus]));
                #plt.title("(PostSynCell: " + str(i_post) + ") Distribution of the delay of spikes after spikes of PreSynCell " + str(i))
                plt.xlabel("delay [ms]");
                plt.gca().invert_xaxis()
                #plt.ylabel("index of neuron in layer 0");
                plt.yticks(range(len(preList)), preList)
            
            
            #tracing poly
            poly_indexs = []
            poly_delays = []
            
            postMax = polyTableTmpConnected.max();
            postArgMax = np.where(polyTableTmpConnected==postMax);
            
            
            if postMax>=polyChainDetectTh:
                for max_in in range(len(postArgMax[0])):
                    poly_indexs.append(preList[postArgMax[0][max_in]]);
                    poly_delays.append(int(postArgMax[1][max_in]));
                    traceDelay(poly_indexs,poly_delays,preList[postArgMax[0][max_in]],int(postArgMax[1][max_in]))
                    plt.subplot(int(subplotDim)+1,2,2);
                    plt.plot(poly_delays,poly_indexs,'*');
                    plt.xlim([max_delay*-0.01,max_delay*1.01]);
                    plt.ylim([-1,nCells]);
                    plt.gca().invert_yaxis()
                    plt.gca().invert_xaxis()
                    #plt.show()
             
            fig_poly.savefig(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/poly"+str(t_min)+"-"+str(t_max)+"/poly_"+str(t_min)+"-"+str(t_max)+"_"+str(i_post)+"_polyLen"+str(len(poly_indexs))+".png");

            
    
    #print "stop";
    #plt.plot(prevSpikeTime[:][1::]);
    
    #plt.show()
    

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
