import numpy as np
import pylab as plt
import pickle
import os;


plotAllSingleCellInfo = 0;

ExperimentName = "practice";
layer = 2; #0: gabor, 1: 1st layer ..

phases = ['FR_0_blank.pkl', 'FR_1_trained.pkl'];

index = 1;
for phase in phases:
    FR=pickle.load(open(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+ExperimentName+"/"+ phase, "rb"));
    numObj = FR.shape[0];
    numTrans = FR.shape[1];
    numLayers = FR.shape[2]-1;
    numCells = FR.shape[3];
    
    FR_layer = FR[:,:,layer,:,:];
    
    FR_norm = (FR-FR[:,:,layer,:,:].min())/(FR[:,:,layer,:,:].max()-FR[:,:,layer,:,:].min());
    FR = FR_norm;
    
    #settings
    numBins =  3;##numTrans;   #can be adjusted
    weightedAnalysis = 0;#exclude the selectivity by not responding to a particular stimulus
    
    #params for multiple cell info (to be implemented)
    multi_cell_analysis = 0; #1 to run multi-cell info analysis
    num_samples = 10;   #num samples from Max-Cells of each stimulus for multi-cell info. analysis
    sampleMulti =1 #4*4;    #default is 1; this is to increase the sampling size to deal with larger layers.
    nc_max = 15;#num_samples*numObj;   # max ensemble size for multi-cell info. analysis
    
    sumPerBin = np.zeros((numCells,numCells,numBins));
    sumPerObj = numTrans;
    sumPerCell = numTrans*numObj;
    IRs = np.zeros((numObj,numCells,numCells));#I(R,s) single cell information
    IRs_weighted = np.zeros((numObj,numCells,numCells));#I(R,s) single cell information
    pq_r = np.zeros(numObj)#prob(s') temporally used in the decoing process
    Ps = 1/numObj   #Prob(s) 
    
    Iss2_Q_avg = np.zeros(nc_max);#average quantized info for n cells; I(s,s')
    Iss2_P_avg = np.zeros(nc_max);## average smoothed info for n cells; I(s,s')
    
    
    print("**Loading data**")
    binMatrix = np.zeros((numCells, numCells, numObj, numBins));# #number of times when fr is classified into a specific bin within a specific objs's transformations
    binMatrixTrans = np.zeros((numCells, numCells, numObj, numBins, numTrans));  #TF table to show if a certain cell is classified into a certain bin at a certain transformation
    for obj in range(numObj):
        print str(obj) + '/' + str(numObj);
        for trans in range(numTrans):
            for row in range(numCells):
                for col in range(numCells):
                    bin = np.around(FR[obj,trans,layer,col,row]*(numBins-1));
                    binMatrix[row,col,obj,bin]=binMatrix[row,col,obj,bin]+1;
                    binMatrixTrans[row,col,obj,bin,trans]=1;
    
    #print binMatrix;
    
    print "** single-cell information analysis **";
    # Loop through all cells to calculate single cell information
    for row in range(numCells):
        for col in range(numCells):
            # For each cell, count the number of transforms per bin
            for bin in range(numBins):
                for obj in range(numObj):
                    sumPerBin[row,col,bin]+=binMatrix[row,col,obj,bin];
    
            # Calculate the information for cell_x cell_y per stimulus
            for obj in range(numObj):
                for bin in range(numBins):
                    Pr = sumPerBin[row,col,bin]/sumPerCell;
                    Prs = binMatrix[row,col,obj,bin]/sumPerObj;
                    if(Pr!=0 and Prs!=0 and Pr<Prs):
                        IRs[obj,row,col]+=(Prs*(np.log2(Prs/Pr)));#*((bin-1)/(numBins-1)); #could be added to weight the degree of firing rates.
                        #IRs(row,col,obj)=IRs(row,col,obj)+(Prs*(log2(Prs/Pr)))*((bin-1)/(numBins-1)); #could be added to weight the degree of firing rates.
                        IRs_weighted[obj,row,col]+=(Prs*(np.log2(Prs/Pr)))*((bin)/(numBins-1)); #could be added to weight the degree of firing rates.
     
    if (weightedAnalysis==1):
        IRs = IRs_weighted;
    
    IRs_reshaped = np.reshape(IRs,(numObj,numCells*numCells));
    IRs_sorted = np.transpose(np.sort(IRs_reshaped));
    reversed_arr = IRs_sorted[::-1]

    plt.subplot(1,2,index)
    plt.plot(reversed_arr);
    plt.ylim([-0.05, np.log2(numObj)+0.05]);
    plt.xlim([0, numCells*numCells])
    index+=1;
    
    
plt.show();
    
