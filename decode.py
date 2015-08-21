def Decode():
    
    for layer in range(0,nLayers):
        
        spikes = spkdetLayers[layer].spike_trains();
        
        spkcount = [len(spikes)][nStim*nTrans]
        spkStimcount = [len(spikes)][nStim]
        
        for i in range(0, len(spikes)-1):
            
            for j in range(0, nStim*nTrans-1):
                
                if (len(spikes[i])==0):
                    condition =   spikes[i]>simulationTime*j && spikes[i]<simulationTime*(j+1)
                else:
                    condition = spikes[i]>simulationTime*j*ms && spikes[i]<simulationTime*(j+1)*ms
                
                spkcount[i][j] = len(numpy.extract(condition,spikes[i]));
                
    



    return 0