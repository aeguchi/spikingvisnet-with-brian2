from Parameters import *


net = Network(collect())

#Creating Gabor input layer
layerG = [];
for theta in range(0,len(thetaList)):
    #layerG.append(NeuronGroup(layerGDim*layerGDim, eqs, threshold='v>1', reset='v = 0'))
    layerG.append(PoissonGroup(layerGDim*layerGDim,numpy.random.rand(layerGDim*layerGDim)*Hz ))
net.add(layerG);

#Creating ExcitLayers
layers = []
for layer in range(0,nLayers):
    #layers.append(NeuronGroup(layerDim*layerDim, eqs, threshold='v>1', reset='v = 0'))
    layers.append(NeuronGroup(layerDim*layerDim, eqn_membran, threshold='v>Vt', reset='v = Vr', refractory=2*ms))
    layers[layer].v = 'Vr'# + rand() * (Vt - Vr)'
    layers[layer].ve = 0*mV
    layers[layer].vi = 0*mV

net.add(layers);

#Creating InhibLayers
inhibLayers = []
for layer in range(0,nLayers):
    #inhibLayers.append(NeuronGroup(inhibLayerDim*inhibLayerDim, eqs, threshold='v>1', reset='v = 0'))
    inhibLayers.append(NeuronGroup(layerDim*layerDim, eqn_membran, threshold='v>Vt', reset='v = Vr', refractory=2*ms))
    inhibLayers[layer].v = 'Vr'# + rand() * (Vt - Vr)'
    inhibLayers[layer].ve = 0*mV
    inhibLayers[layer].vi = 0*mV
net.add(inhibLayers);


#Connecting neurons in Gabor input layer and neurons in the first ExcitLayer
connGtoInput = []
for theta in range(0,len(thetaList)):
    
    connGtoInput.append(Synapses(layerG[theta], layers[0],'plastic : 1 (shared)', pre='ve += 3*we* plastic'));
    connGtoInput[theta].connect(True, p=0.2)
#connGtoInput[theta].connect(sqrt((i%layerGDim - j%layerGDim)**2 + (i/layerGDim - j/layerGDim)**2) < 0.1 * layerGDim)

#connGtoInput.append(Synapses(layerG[theta], layers[0], 'w:1', pre='ve += w*we'));
#connGtoInput[theta].w = 1;
#connGtoInput.append(Synapses(layerG[theta], layers[0], pre='ve += 10*we'));
#connGtoInput[theta].connect(True, p=0.1)
net.add(connGtoInput)


#Connecting neurons between layers
connBottomUp = []
connTopDown = []
for layer in range(0,nLayers-1):
    #connBottomUp.append(Synapses(layers[layer],layers[layer+1], 'w:1',pre='ve += 10*we'))
    connBottomUp.append(Synapses(layers[layer],layers[layer+1], eqs_stdpSyn, eqs_stdpPre ,eqs_stdpPost))
    connBottomUp[layer].connect(True, p=0.1)
    connBottomUp[layer].w[:,:] = rand()*Apre
    connBottomUp[layer].delay[:,:] = rand()*10*ms
    
    connTopDown.append(Synapses(layers[layer+1],layers[layer], eqs_stdpSyn, eqs_stdpPre ,eqs_stdpPost))
    #connTopDown.append(Synapses(layers[layer+1],layers[layer], 'w:1',pre='ve += 10*we'))
    connTopDown[layer].connect(True, p=0.1);
    connTopDown[layer].w[:,:] = rand()*Apre
    connTopDown[layer].delay[:,:] = rand()*10*ms

net.add(connBottomUp)
net.add(connTopDown)

#Connecting neurons within layers
connExIn = []
connInEx = []
connRecIn = []
connRecEx = []
for layer in range(0,nLayers):
    connExIn.append(Synapses(layers[layer],inhibLayers[layer], '''w:1 
                                                                plastic : 1 (shared)''',pre='ve += we* plastic'))
    #connExIn[layer].w = 1;
    connExIn[layer].connect(True, p=0.2)
    connExIn[layer].delay[:,:] = rand()*2*ms
    
    connInEx.append(Synapses(inhibLayers[layer], layers[layer],'''w:1
                                                                plastic : 1 (shared)''', pre='vi += 3*wi* plastic'))
    #connInEx[layer].w = 1;
    connInEx[layer].connect(True, p=0.2)
    connInEx[layer].delay[:,:] = rand()*2*ms
net.add(connExIn)
net.add(connInEx)



spikesG=[]
for theta in range(0,len(thetaList)):
    spikesG.append(SpikeMonitor(layerG[theta]))
net.add(spikesG)

# tmp = layerG[0];
# testSpikes = SpikeMonitor(tmp);

spkdetLayers = []
for layer in range(0,nLayers):
    spkdetLayers.append(SpikeMonitor(layers[layer]))
#spkdetLayers.append(StateMonitor(layers[layer], 'v', record=True))
net.add(spkdetLayers)

spkdetInhibLayers = []
for layer in range(0,nLayers):
    spkdetInhibLayers.append(SpikeMonitor(inhibLayers[layer]))
net.add(spkdetInhibLayers)


def plasticON(on):
    if (on==1):
        for theta in range(0,len(thetaList)):
            connGtoInput[theta].plastic = True

        for layer in range(0,nLayers):
            connExIn[layer].plastic= True
            connInEx[layer].plastic= True

        for layer in range(0,nLayers-1):
            connBottomUp[layer].plastic= True
            connTopDown[layer].plastic= True

    elif (on==0):
        for theta in range(0,len(thetaList)):
            connGtoInput[theta].plastic = False
        
        for layer in range(0,nLayers):
            connExIn[layer].plastic= False
            connInEx[layer].plastic= False
        
        for layer in range(0,nLayers-1):
            connBottomUp[layer].plastic= False
            connTopDown[layer].plastic= False
    return 0
