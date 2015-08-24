from itertools import groupby

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
        
    spikes = [nStim*nTrans][nLayers][LayerDim][2]
    tmp = [trialNb*TestNb]
                        
        for k in range(nStim*nTrans):
            for l in range(nLayers):
                for n in range(LayerDim*LayerDim):
                    for trial in range(trialNb) :
                        for test in range(testNb):

                            if (len(FR_t[l][trial][test][n])==0):
                                condition = FR_t[l][trial][test][n]>simulationTime*k && FR_t[l][trial][test][n]<simulationTime*(k+1)
                            else:
                                condition = FR_t[l][trial][test][n]>simulationTime*k*ms && FR_t[l][trial][test][n]<simulationTime*(k+1)*ms

                            tmp[trial*trialNb + test*testNb] = len(numpy.extract(condition,FR_t[layer][trial][test][n]))
                
                    tmp=numpy.sort(tmp)
                    Pr=[len(list(group)) for key, group in groupby(tmp)]
                    spikes[k][l][n][1] = tmp
                    spikes[k][l][n][2] =Pr


    
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