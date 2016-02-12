from brian2 import ms;
import numpy as np;
import pylab as plt;
import pickle
import os;
import errno

experimentName = 'test_V1_norm1_G2L3_I2E_in';
SpikeData = ['0_spikes_e.pkl', '0_spikes_i.pkl', '0_pikes_b.pkl', '0_pikes_g.pkl'];
spikes_e = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/" + SpikeData[0], "rb"));
netState_L2L = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/0_netStates_L2L.pkl", "rb"));

preSynConn = netState_L2L[0]['i_pre'];
postSynConn = netState_L2L[0]['i_post'];
weightList = netState_L2L[0]['w'];
delayList = netState_L2L[0]['delay'];




# analysing cell index==0 for test development
nCells = 100;
t_min = 0;#319000;
t_max = 5000;
max_delay = 20;
gmax = 2.0

try:
    os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/hist"+str(t_min)+"-"+str(t_max));
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise

fig_hist = plt.figure(1 , figsize=(40, 40),dpi=500);
for i_post in range(100):
    print i_post
    plt.clf();

    prevSpikeTime = np.zeros([nCells, len(spikes_e[1][i_post])]);  # pre cell index, number of spikes
    spikeTime_post = spikes_e[1][i_post];
    cond_post1 = spikeTime_post < t_max;
    cond_post2 = t_min < spikeTime_post;
    cond_post = cond_post1 & cond_post2;
    spikeTime_post_ext=np.extract(cond_post, spikeTime_post)

    if(len(spikeTime_post_ext)>0):
        count = 0;
        for t_post in spikeTime_post_ext:
            for i_pre in range(100):
                if len(spikes_e[0][i_pre]) > 0:
                    # print i_pre
                    spikeTime_pre = spikes_e[0][i_pre];
                    spikeTime_pre = spikeTime_pre[::-1];
                    cond_pre1 = t_post-max_delay < spikeTime_pre;
                    cond_pre2 = spikeTime_pre < t_post;
                    cond = cond_pre1 & cond_pre2;
                    spikeTime_pre_ext = np.extract(cond, spikeTime_pre);
                    if (len(spikeTime_pre_ext) > 0):
                        t_back = t_post - spikeTime_pre_ext[0];
                        prevSpikeTime[i_pre][count] = t_back;
            count += 1;
            #print count
        
        for i_pre in range(100):
            plt.subplot(10,10,i_pre);
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

            tmp = np.extract(prevSpikeTime[i_pre]!=0,prevSpikeTime[i_pre]);
            if len(tmp)>1:
                plt.hist(tmp,bins=10,range = (0,max_delay),color=c);
                
            if len(preList)>0:
                plt.hold(True);
                plt.plot(delay[0],len(spikeTime_post_ext)/2,'*');
                
            plt.ylim(len(spikeTime_post_ext)*-0.1,len(spikeTime_post_ext)*1.1)
            plt.xlim(-0.01*max_delay,max_delay*1.01);
            
            
        fig_hist.savefig(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/hist"+str(t_min)+"-"+str(t_max)+"/hist_"+str(t_min)+"-"+str(t_max)+"_"+str(i_post)+".png");

    
    
    
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
