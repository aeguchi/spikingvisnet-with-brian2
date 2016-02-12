from brian2 import ms;
import numpy as np;
import pylab as plt;
import pickle
import os;
import errno

experimentName = 'test_V1_norm1_itr500';
#SpikeData = ['0_spikes_e.pkl', '0_spikes_i.pkl', '0_pikes_b.pkl', '0_pikes_g.pkl'];
SpikeData = ['501_spikes_e.pkl', '501_spikes_i.pkl', '501_pikes_b.pkl', '501_pikes_g.pkl'];
spikes_e = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/" + SpikeData[0], "rb"));
#netState_L2L = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/0_netStates_L2L.pkl", "rb"));
netState_L2L = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/502_netStates_L2L.pkl", "rb"));

preSynConn = netState_L2L[0]['i_pre'];
postSynConn = netState_L2L[0]['i_post'];
weightList = netState_L2L[0]['w'];
delayList = netState_L2L[0]['delay'];

polyAnalysisOn = True
polyHist = False

nBins = 10;
# analysing cell index==0 for test development
nCells = 100;
t_min = 320000;#319000;
t_max = 325000;
max_delay = 20;
gmax = 2.0

try:
    os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/hist"+str(t_min)+"-"+str(t_max));
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise

#fig_hist = plt.figure(1 , figsize=(40, 40),dpi=500);
#fig_poly = plt.figure(2 , figsize=(40, 40),dpi=100);

for i_post in range(100):
    print i_post
    #plt.clf();

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
            polyTable = np.zeros([nCells,nCells,nCells]);#post,pre,diff
            for i in range(count):#for each post synaptic spike
                for cell_in in range(nCells):            
                    tmp = prevSpikeTime[:,i];  #store delay of each presynaptic spike (before the ith post synaptic spike)
                    if tmp[cell_in]<max_delay+1:
                        sorted = np.sort(tmp)       #sorted the delay of each input cell (before the ith post synaptic spike)
                        sorted_i = np.argsort(tmp)  #sorted index
                        cond_sorted1 = sorted<max_delay+1; #eliminated all the stored delay bigger than max_delay
                        cond_sorted2 = tmp[cell_in]<sorted
                        cond_sorted = cond_sorted1&cond_sorted2;
                        ext_sorted = np.extract(cond_sorted,sorted);
                        ext_sorted_i = np.extract(cond_sorted,sorted_i);
                        t_backs=ext_sorted-tmp[cell_in];
                        #print str(cell_in) + " " + str(t_backs);
                        for order_in in range(len(ext_sorted)):
                            t_back = int(t_backs[order_in]/1);
                            polyTable[cell_in][ext_sorted_i[order_in]][t_back] += 1;        
            
            for i in range(nCells):
                fig_poly = plt.figure(2 , figsize=(40, 40),dpi=100);
                #print polyTable[i]
                plt.clf()
                plt.imshow(polyTable[i],interpolation='none');
                plt.colorbar()
                plt.title("(PostSynCell: " + str(i_post) + ") Distribution of the delay of spikes after spikes of PreSynCell " + str(i))
                plt.xlabel("delay [ms]");
                plt.ylabel("index of neuron in layer 0");
                plt.show(fig_poly)    
    
    
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
