using Pkg
dir = pwd()
Pkg.activate("$(dir)/.")
using BioNeuralMass, Sundials, DifferentialEquations, Plots, Waveforms

nm = BioNM() # Sets up neural mass struct

# Set up cortical neural area
hp = HyperParam(synapseoff = false, ratio = 0.8)
pop1 = NeuralPopSoma{Cortex, Excitatory}() # Excitatory cortical population
pop2 = NeuralPopSoma{Cortex, Inhibitory}() # Inhibitory cortical population 
cortex_ = NeuralArea(hp,pop1,pop2); # Cortex region
Par(hp,cortex_) # computes impermeants and leak conductances

# Set up Thalamic neural area
pop1_T = NeuralPopSoma{Thalamus, Excitatory}() # Excitatory thalamic population
pop2_T = NeuralPopSoma{Thalamus,Inhibitory}() # Inhibitory thalamic population
thalamus_ = NeuralArea(hp,pop1_T,pop2_T) # Thalamus region
Par(hp,thalamus_) # computes impermeants and leak conductances

# Set up neural mass
# nm.conn is a matrix of connections, (i,j) refers to strength of current (j-> i).
nm.conn = [ 0.0 2  0.5 0.0 ; 0.3 0.0 0.1 0.0; 0.3 0.0  0.3 5; 0.2 0 0.5 2.5] #Thalamus first
nm.conn = [ 0.0 2  0.5 0.0 ; 0.3 0.0 0.1 0.0; 0.5 0.0  0.3 5; 0.2 0 0.5 2.5] #Thalamus first


# Thalamus baseline conditions 
thalamus_.pop1.syn_th = 0.2 
thalamus_.pop1.syn_act = 1.25
thalamus_.pop1.syn_deact = 0.3
thalamus_.pop2.syn_th = 0.5
thalamus_.pop2.syn_act = 0.5
thalamus_.pop2.syn_deact = 0.003

# Moderate ischemia weakens cortical activation, makes it identical to cortex
fac = 1
cortex_.pop1.syn_act  = fac*cortex_.pop1.syn_act   
cortex_.pop1.syn_deact= fac*cortex_.pop1.syn_deact
cortex_.pop2.syn_act  = fac*cortex_.pop2.syn_act  
cortex_.pop2.syn_deact= fac*cortex_.pop2.syn_deact

# Hyperparameterd
hp.excite = [20, 1000,5*60*1e3,1e8] # [Current strength, start time, end time, large number]
# hp.excite = missing # no stimulation
hp.perc = 1# Available energy ∈ [0,1], 0 -> complete ED, 1 -> No ED
hp.tstart = 1*60*1e3 # ED start time
hp.tend =  3*60*1e3 # ED end time
hp.tfinal = 1e4 # Simulation end time
hp.beta1 = 5 # ED onset rate (higher is faster)
hp.beta2 = 5 # ED offset rate (higher is faster)
hp.saveat = 0.1 # save every at every x milliseconds
hp.O2e_th_vATP = 1.5 # vATPase oxygen threshold
hp.O2e_th_NKA = 1.5 # NKA oxygen threshold

nm.areas = [thalamus_,cortex_] # Order of areas matters! First area is the one stimulated. Take care of nm.conn!
nm.hp = hp


#fac2 = range(0,2,length=50)
#fac1 = range(0.2,1.2,length=25)
#fac = range(0,2,length=20)
#conn41 = range(0.0,0.04,length=100)
#nm.conn[3,1] = 0.45
#nm.conn[4,1] = 0.18
# nm.conn[4,1] = 0.1255 

# solve(nm,saveat=hp.saveat,reltol=1e-9,abstol=1e-9) # Solve system using solve(nm)
# solve(nm,CVODE_BDF(),saveat=hp.saveat,reltol=1e-7,abstol=1e-7) # Solve system using solve(nm)
# plot_syn(nm,1,"time (ms.)","ModerateED") # Plot synaptic currents, change second argument to 60*1000 to plot in min.
# plot(nm,"Baseline") # Plots EEG and ion dynamics for both regions and saves as two files.
