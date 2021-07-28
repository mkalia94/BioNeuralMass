using Pkg
dir = pwd()
Pkg.activate("$(dir)/.")
using BioNeuralMass, Sundials, DifferentialEquations, Plots, Waveforms

nm = BioNM()

# Set up cortical neural area
hp = HyperParam(synapseoff = false, ratio = 0.8)
pop1 = NeuralPopSoma{Cortex, Excitatory}()
pop2 = NeuralPopSoma{Cortex, Inhibitory}()
cortex_ = NeuralArea(hp,pop1,pop2);
Par(hp,cortex_)

# Set up Thalamic neural area
pop1_T = NeuralPopSoma{Thalamus, Excitatory}()
pop2_T = NeuralPopSoma{Thalamus,Inhibitory}()
thalamus_ = NeuralArea(hp,pop1_T,pop2_T)
Par(hp,thalamus_)

thalamus_.pop1.syn_act *= 1e-1
thalamus_.pop2.syn_act *= 1e-1

thalamus_.pop1.syn_th = 0.5
thalamus_.pop2.syn_th = 0.9


# Set up neural mass
nm.conn = [0.0 1 2.0 0.0; 0.5 0.0 1.5 0.0; 3.0 0.0 2.0  3.0; 3.0 0.0 1.0 0.5] #Thalamus first 
nm.conn = [2.0 3.0  0.01 0.0; 1.0 0.5 0.01 0.0; 0.01 0.0 0.0 1; 0.01 0.0 0.5 0.0] #Cortex first 

hp.excite = [20.0, 1000.0,1000.0 + 2*60*1000, 1e10]
hp.perc = 0.2
hp.tstart = 1000
hp.tend = 1000 + 2*60*1000
hp.tfinal = 3*60*1000
hp.beta1 = 1
hp.beta2 = 1
hp.saveat = 1

solve(nm,hp,[cortex_, thalamus_],saveat=hp.saveat,reltol=1e-9,abstol=1e-9)
syn_curr = nm([cortex_,thalamus_],hp,"syn_curr")
plot(syn_curr,layout = (2,2))
