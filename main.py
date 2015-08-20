
from imageImport import *

from PlotSim import *
import sys
import  glob
#import nest
import os
#from nest import raster_plot
#import nest.topology as topp




#http://brian2.readthedocs.org/en/latest/resources/tutorials/1-intro-to-brian-neurons.html




### Training ###
#trainingImages = sorted(glob.iglob("/images/training/*.png"))
trainingImages = sorted(glob.iglob(os.path.split(os.path.realpath(__file__))[0] + "/images/training/*.png"))

store('initialized')    #randomised network

if not os.path.exists(os.path.split(os.path.realpath(__file__))[0] + "/randomised/"):
    os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/randomised/")
if not os.path.exists(os.path.split(os.path.realpath(__file__))[0] + "/trained/"):
    os.makedirs(os.path.split(os.path.realpath(__file__))[0] + "/trained/")

img_fns = []
for img_fn in trainingImages:
    img_fns.append(img_fn)

if len(img_fns)!=nStim*nTrans:
    print 'Error: the number of images files does not match',len(img_fns);
    sys.exit(1)

index_img=0;
for img_fn in img_fns:
    #print img_fn
    #load gabor filtered ImageSurface
    img = cv2.imread(img_fn)
    if img is None:
        print 'Failed to load image file:', img_fn
        sys.exit(1)
    
    filters = build_filters()
    res = process(img, filters)
    
    #convert from gabor filtered inputs to spikes
    #normalize
    res_norm=res/numpy.max(res);
    res_norm=1-res_norm;
    
    for index_filter in range(0,len(thetaList)):
        r = numpy.reshape(numpy.mean(res_norm[index_filter],axis=2),(layerGDim*layerGDim));
        #print r
        layerG[index_filter].rates= r * Rmax;    #To be fixed
        #print layerG[index_filter].rates
        
    
    
    
    net.run(simulationTime*ms)
    plotSim(img,res_norm,res,index_img)

    index_img+=1;
    
#To-DO: trainNetworkWith
if (plotLayer==1):
#         plt.figure(2);    
#         plt.subplot(nLayers+1,3,1);
#         plt.imshow(img,interpolation='none');
#         plt.title('Input')
    
    for layer in range(0,nLayers):
        #headNodeIndex = layers[layer][0]

        #plot spike raster
        ax=plt.subplot(nLayers,4,(nLayers-layer-1)*4+1);
        
        plt.title(layer)
        raster_plot(spkdetLayers[layer])
        pylab.ylim([0,layerGDim*layerGDim-1])
        pylab.xlim([simulationTime*index_img,simulationTime*index_img+1])
        #ax.get_yaxis().set_visible(False)
        #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
        #ax.set_ylim([headNodeIndex+1, headNodeIndex+(layer1Dim*layer1Dim)]);
        
        
        #dSD =nest.GetStatus(spkdetInhibLayers[layer],keys='events')[0];
            
        #plot spike raster
        ax=plt.subplot(nLayers,4,(nLayers-layer-1)*4+3);
        plt.title(layer)
        raster_plot(spkdetInhibLayers[layer])
        #ax.get_yaxis().set_visible(False)
        #ax.set_xlim([(index_img)*simulationTime, (index_img+1)*simulationTime])
        #ax.set_ylim([headNodeIndex+1, headNodeIndex+(inhibLayer1Dim*inhibLayer1Dim)]);
        
    plt.show()
    
#Testing
