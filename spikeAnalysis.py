from brian2 import ms;
import numpy as np;
import pylab as plt;
import pickle
import os;

experimentName = 'poly4';
SpikeData = ['Spikes_e.pkl', 'Spikes_i.pkl','Spikes_b.pkl','Spikes_g.pkl'];
tmp=pickle.load(open(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/"+ SpikeData[0], "rb"));

#spikes[layer][cellIndex]
layer = 1;
spikes = tmp[layer];
endingTime = 7200; 
print endingTime
plt.subplot(3,1,1);
plt.hold(True)
for i in range(100):
    plt.plot(spikes[i]/ms,np.ones(len(spikes[i]))*i,'.');
plt.xlim([endingTime-3000,endingTime-2000])

plt.subplot(3,1,2);
plt.hold(True)
for i in range(100):
    plt.plot(spikes[i]/ms,np.ones(len(spikes[i]))*i,'.');
plt.xlim([endingTime-2000,endingTime-1000])

plt.subplot(3,1,3);
plt.hold(True)
for i in range(100):
    plt.plot(spikes[i]/ms,np.ones(len(spikes[i]))*i,'.');
plt.xlim([endingTime-1000,endingTime])

plt.show()