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

# Set up neural mass
nm.conn = [0.0 1 2.0 0.0; 0.5 0.0 1.5 0.0; 3.0 0.0 2.0  3.0; 3.0 0.0 1.0 0.5] #Thalamus first 
nm.conn = [ 0.0 2  0.5 0.0 ; 0.3 0.0 0.1 0.0; 0.3 0.0  0.3 5; 0.2 0 0.5 2.5] #Cortex first
# nm.conn = [ 0.7 5  0.0 0.0 ; 0.5 2.5 0.0 0.0; 0.5 0.0  0.7 5; 0.2 0 0.5 2.5] #Cortex first

# Thalamic edits
thalamus_.pop1.syn_th = 0.2
thalamus_.pop1.syn_act = 1.25
thalamus_.pop1.syn_deact = 0.3
thalamus_.pop2.syn_th = 0.5
thalamus_.pop2.syn_act = 0.5
thalamus_.pop2.syn_deact = 0.003


hp.excite = [20, 60*1000.0,60*1000 + 2*60*1000, 1e10]
hp.excite = missing
hp.perc = 0.7
hp.tstart = 60*1000
hp.tend = 60*1000 + 2*60*1000
hp.tfinal = 4*60*1000
hp.beta1 = 1
hp.beta2 = 1
hp.saveat = 0.5


nm.areas = [thalamus_,cortex_]
nm.hp = hp

solve(nm,saveat=hp.saveat,reltol=1e-9,abstol=1e-9)
plot_syn(nm,60*1000,"time (min.)","Perc10")
