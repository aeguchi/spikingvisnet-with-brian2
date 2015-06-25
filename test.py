import nest
import pylab
from random import random

nest.Models()

# neuron = nest.Create("iaf_neuron")
# 
# nest.GetStatus(neuron)
# 
# #nest.SetStatus(neuron, {"I_e": 376.0})
# 
# multimeter = nest.Create("multimeter")
# nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})
# 
# 
# spikedetector = nest.Create("spike_detector",
# params={"withgid": True, "withtime": True})
# 
# nest.Connect(multimeter, neuron)
# nest.Connect(neuron, spikedetector)
# 
# 
# 
# 






preNeurons = nest.Create("iaf_neuron",10)
#nest.SetStatus(preNeurons, {"I_e": 376.0})
postNeurons = nest.Create("iaf_neuron",10)

multimeter = nest.Create("multimeter",10)
nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})

# spikedetector = nest.Create("spike_detector", 10)
# nest.SetStatus(spikedetector, {"withgid": True, "withtime": True})
spikedetector = nest.Create("spike_detector", 10, params={"withgid": True, "withtime": True})


#noise generator
noiseEx = nest.Create("poisson_generator", 10)
nest.SetStatus(noiseEx, [{"rate": 80000.0*random()},{"rate": 80000.0*random()},{"rate": 80000.0*random()},{"rate": 80000.0*random()},{"rate": 80000.0*random()},{"rate": 80000.0*random()},{"rate": 80000.0*random()},{"rate": 80000.0*random()},{"rate": 80000.0*random()},{"rate": 80000.0*random()}])

syn_dict_ex = {"weight": 1.2}

noiseIn = nest.Create("poisson_generator", 10)
nest.SetStatus(noiseIn, [{"rate": 15000.0*random()},{"rate": 15000.0*random()},{"rate": 15000.0*random()},{"rate": 15000.0*random()},{"rate": 15000.0*random()},{"rate": 15000.0*random()},{"rate": 15000.0*random()},{"rate": 15000.0*random()},{"rate": 15000.0*random()},{"rate": 15000.0*random()}])
syn_dict_in = {"weight": -0.5}





#establishing connections
nest.Connect(preNeurons, postNeurons, syn_spec = {"weight":20.0})
#nest.Connect(preNeurons, postNeurons, "all_to_all", {"weight":1.0})#fully connected
nest.Connect(noiseIn, preNeurons, syn_spec=syn_dict_in)
nest.Connect(noiseEx, preNeurons, syn_spec=syn_dict_ex)
nest.Connect(multimeter, postNeurons)
nest.Connect(postNeurons,spikedetector)





nest.Simulate(1000.0)

dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["V_m"]
ts = dmm["events"]["times"]

pylab.figure(1)
pylab.plot(ts, Vms)


dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
pylab.figure(2)
pylab.plot(ts, evs, ".")







pylab.show()


var = raw_input();