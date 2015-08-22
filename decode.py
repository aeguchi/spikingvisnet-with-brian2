def Decode():
    
    FR_t=[nlayers][trialNb][testNb]
    FR_b=[nlayers][trialNb]
    
    
    for layer in range(0,nLayers):
        
        for trial in range(trialNb) :
            
            filep= os.path.split(os.path.realpath(__file__))[0] + '/FR/' + 'FR_b_Ex'+str(layer)+'_trial'+srt(trial)
            with open(filep,'w') as FR
                FR_b[layer][trial] = json.load(FR)
            
            for test in range(testNb):
                filep= os.path.split(os.path.realpath(__file__))[0] + '/FR/' + 'FR_b_Ex'+str(layer)+'_trial'+srt(trial)+'_test'+str(test)
                with open(filep,'w') as FR
                    FR_t[layer][trial][test] = json.load(FR)
        
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