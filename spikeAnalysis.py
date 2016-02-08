from brian2 import ms;
import numpy as np;
import pylab as plt;
import pickle
import os;

experimentName = 'poly1000';
SpikeData = ['Spikes_e.pkl', 'Spikes_i.pkl', 'Spikes_b.pkl', 'Spikes_g.pkl'];
spikes_e = pickle.load(open(os.path.split(os.path.realpath(__file__))[0] + "/Results/" + experimentName + "/" + SpikeData[0], "rb"));


# analysing cell index==0 for test development
nCells = 100;
t_max = 1000;
i_post = 1;
max_delay = 20;

fig_hist = plt.figure(0 , figsize=(20, 20),dpi=500);
for i_post in range(100):
    print i_post
    plt.clf();

    prevSpikeTime = np.zeros([nCells, len(spikes_e[1][i_post])]);  # pre cell index, number of spikes
    spikeTime_post = np.sort(spikes_e[1][i_post] / ms);
    
    count = 0;
    for t_post in (np.extract(spikeTime_post < t_max, spikeTime_post)):
        for i_pre in range(100):
            if len(spikes_e[0][i_pre]) > 0:
                # print i_pre
                spikeTime_pre = np.sort(spikes_e[0][i_pre] / ms);
                spikeTime_pre = spikeTime_pre[::-1];
                cond1 = t_post-max_delay < spikeTime_pre;
                cond2 = spikeTime_pre < t_post;
                cond = cond1 & cond2;
                extractedSpikeTime = np.extract(cond, spikeTime_pre);
                if (len(extractedSpikeTime) > 0):
                    t_back = t_post - extractedSpikeTime[0];
                    prevSpikeTime[i_pre][count] = t_back;
        count += 1;
        #print count
    
    for i_pre in range(100):
        plt.subplot(10,10,i_pre);
        tmp = np.extract(prevSpikeTime[i_pre]!=0,prevSpikeTime[i_pre]);
        if len(tmp)>1:
            plt.hist(tmp,range = (0,20));
        plt.ylim(-1,11)
        plt.xlim(-1,21);
    fig_hist.savefig(os.path.split(os.path.realpath(__file__))[0] + "/Results/"+experimentName+"/hist"+str(i_post)+".png");

    
    
    
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
