import numpy as np
import pylab as plt
import pickle
import os;
#from Parameters import *


class InfoAnalysis(object):
    
    def __init__(self,borrowed_globals):
        globals().update(borrowed_globals);

    def singleCellInfoAnalysis(self):
        plotAllSingleCellInfo = 0;
        
        #ExperimentName = "BO_single";
        #layer = 2; #0: bindingLayer, 1: 1st layer ..
        plotMode = 1; #0:normal 1:max
        
        phases = ['FR_0_blank.pkl', 'FR_1_trained.pkl'];
        linestyles = ['-', '--'];#, '-.', ':']
        
        
        
        FR=pickle.load(open(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/"+ phases[0], "rb"));
        numCells = FR.shape[3];
        
        infos = np.zeros((2,numCells*numCells));
        
        performanceMeasure = 0.0;
        
        for layer in range(3):
            index = 1;
            for phase in phases:
                FR=pickle.load(open(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/"+ phase, "rb"));
                numObj = FR.shape[0];
                numTrans = FR.shape[1];
                numLayers = FR.shape[2]-1;
                numCells = FR.shape[3];
                
                FR_layer = FR[:,:,layer,:,:];
                
                if FR[:,:,layer,:,:].max()>0.001:
                    FR_norm = (FR-FR[:,:,layer,:,:].min())/(FR[:,:,layer,:,:].max()-FR[:,:,layer,:,:].min());
                    FR = FR_norm;
                
                #settings
                numBins =  3;##numTrans;   #can be adjusted
                weightedAnalysis = 1;#exclude the selectivity by not responding to a particular stimulus
                
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
                if plotMode==1:
                    IRs_reshaped = np.max(IRs_reshaped,axis=0);
                
                IRs_sorted = np.transpose(np.sort(IRs_reshaped));
                reversed_arr = IRs_sorted[::-1]
            
            
                infos[index-1] = reversed_arr;
            
                if plotMode==0:
                    plt.subplot(1,2,index)
                    plt.plot(reversed_arr, color='k');
                    plt.ylim([-0.05, np.log2(numObj)+0.05]);
                    plt.xlim([0, numCells*numCells])
                    plt.title(phase);
                    
                index+=1;
                
                
                if (layer == 2 and phase == 'FR_1_trained.pkl'):
                    numInfoCalc = 100;
                    performanceMeasure = -1*np.sum(reversed_arr[0:numInfoCalc]);
                    f = open(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/performance.txt","w");
                    f.write(str(performanceMeasure));
                    f.close();
                    
                    
                    
                
            if plotMode==1:
                plt.subplot(1,3,layer+1)
                plt.plot(np.transpose(infos[0]), linestyle='--', color='k');
                plt.hold(True);   
                plt.plot(np.transpose(infos[1]), linestyle='-', color='k');
                plt.ylim([-0.05, np.log2(numObj)+0.05]);
                plt.xlim([0, numCells*numCells])
                if layer==0:
                    plt.title("Binding layer");
                else:
                    plt.title("Layer " + str(layer));
            
            
        

        plt.savefig(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/singleCellInfo.png");
        plt.savefig(os.path.split(os.path.realpath(__file__))[0] +"/Results/"+experimentName+"/singleCellInfo.eps");
        #plt.show();  
