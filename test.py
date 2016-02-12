from Parameters import *
import brian2 as br
import pylab as plt
import numpy as np

plt.figure(0,figsize=(30, 10),dpi=100)

conductanceConst_I2E = 2.0









numInputCells = 10;
numConn = 1;
#input
layerG=br.PoissonGroup(numInputCells, Rmax)  # rates
# Excitatory layers
#layerG=br.SpikeGeneratorGroup(numInputCells, 't % (1.0 / (br.rand(50)*Hz)) == 0*ms');


# Excitatory layers
layerEx=br.NeuronGroup(1,  # N neurons
                       eqn_membranEx,  # differential equations
                       threshold='v>Vth_ex',  # spike condition
                       reset='''v = V0_ex
                                ge = 0
                                gi = 0''',  # code to execute on reset
                       refractory=refractoryPeriod  # length of refractory period
                       );

layerEx.v = 'V0_ex'  # + br.rand() * (Vt - Vr)'
layerEx.ge = 0
layerEx.gi = 0
#layerEx.Iext = 25*mV;


# Excitatory layers
layerIn=br.NeuronGroup(1,  # N neurons
                       eqn_membranIn,  # differential equations
                       threshold='v>Vth_in',  # spike condition
                       reset='''v = V0_in
                                ge = 0
                                gi = 0''',  # code to execute on reset
                       refractory=refractoryPeriod  # length of refractory period
                       );

layerIn.v = 'V0_in'  # + br.rand() * (Vt - Vr)'
layerIn.ge = 0
layerIn.gi = 0
#layerEx.Iext = 25*mV;

#conn1=br.Synapses(layerG, layerEx,eqs_Syn, pre=eqs_ExPre);

conn1=br.Synapses(layerG, layerEx,eqs_stdpSyn,pre = eqs_stdpPre, post = eqs_stdpPost);
gmax = 2.0;
# conn1.connect('i==0');
# conn1.delay[0] = 0*ms
# 
# conn1.connect('i==1');
# conn1.delay[1] = 0*ms
conn1.connect(True,n=numConn);
#conn1.w[:,:] = 0.1;
conn1.plastic = True;
conn1.w[:,:]=gmax/2;
conn1.delay[:,:]=range(numConn)*ms;




conn2=br.Synapses(layerG,layerIn,eqs_Syn,pre=eqs_InPre);
conn2.connect(True);
conn2.delay[:,:]=10*ms;
conn2.w[:,:]=conductanceConst_I2E;


spkMonLayerG = br.SpikeMonitor(layerG);
spkMonLayerEx = br.SpikeMonitor(layerEx);
sttMonLayerEx = br.StateMonitor(layerEx,'v', record=True,when='before_resets');

#popMonLayerEx = br.PopulationRateMonitor(layerG);

simTime = 100;
nEp = 4
for i in range(nEp):
    print i
    plt.subplot(3,nEp,nEp*0+i+1)
    states1= conn1.get_states(['i_pre','i_post','w']);
    plt.hist(states1['w'],range=[0,gmax]);
    plt.xlim([gmax*-0.01,gmax*1.01]);
    plt.ylim([numInputCells*-0.1,numInputCells*1.1]);
    

    br.run(simTime*ms);

    plt.subplot(3,nEp,nEp*1+i+1);
    plt.plot(sttMonLayerEx.t/ms,sttMonLayerEx.v[0]);
    plt.xlim([simTime*i,simTime*(i+1)])
    plt.ylim([V0_ex*1.01,Vth_ex*0.99])
    
    plt.subplot(3,nEp,nEp*2+i+1);
    plt.hold(True);
    for j in range(numConn):
        plt.plot(spkMonLayerG.t/ms+j,spkMonLayerG.i,'.');
        plt.xlim([simTime*i,simTime*(i+1)])
        plt.ylim([-1,numInputCells*numConn*1.1])
    #plt.hold(True);
    plt.plot(spkMonLayerEx.t / ms,spkMonLayerEx.i+0.5,'*',c='r');
    
#     plt.subplot(4,nEp,nEp*3+i+1);
#     window = 20*ms
#     window_length = int(window/br.defaultclock.dt)
#     cumsum = np.cumsum(np.insert(popMonLayerEx.rate, 0, 0))
#     binned_rate = (cumsum[window_length:] - cumsum[:-window_length]) / window_length
#     plt.bar(popMonLayerEx.t[window_length-1:]/ms, binned_rate);
    
#     plt.xlim([simTime*i+window/ms,simTime*(i+1)])
    
    
plt.show()