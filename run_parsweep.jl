
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
# nm.conn = [0.0 1 2.0 0.0; 0.5 0.0 1.5 0.0; 3.0 0.0 2.0  3.0; 3.0 0.0 1.0 0.5] #Thalamus first 
nm.conn = [ 0.0 2  0.5 0.0 ; 0.3 0.0 0.1 0.0; 0.5 0.0  0.3 5; 0.2 0 0.5 2.5] #Thalamus first
# nm.conn = [ 0.7 5  0.0 0.0 ; 0.5 2.5 0.0 0.0; 0.5 0.0  0.7 5; 0.2 0 0.5 2.5] #Cortex first

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


hp.excite = [20, 1000,20*60*1e3, 1e8] # [Current strength, start time, end time, large number]
# hp.excite = missing # no stimulation
hp.perc = 0.7# Available energy ∈ [0,1], 0 -> complete ED, 1 -> No ED
hp.tstart = 1*60*1e3 # ED start time
hp.tend =  3*60*1e3 # ED end time
hp.tfinal = 5*60*1e3 # Simulation end time
hp.beta1 = 5 # ED onset rate (higher is faster)
hp.beta2 = 5 # ED offset rate (higher is faster)
hp.saveat = 4 # save every at every x milliseconds


nm.areas = [thalamus_,cortex_] # Order of areas matters! First area is the one stimulated. Take care of nm.conn!
nm.hp = hp

fac = [0.3; 0.6; 1]
conn = deepcopy(nm.conn)

for i in 1:length(fac)
    for j in 1:length(fac)
        for k in 1:length(fac)
            for l in 1:length(fac)
                nm.hp.perc = 1
                nm.hp.tfinal = 2*1e4
                nm.hp.saveat = 1
                nm.conn[1,3] = fac[i]*conn[1,3] 
                nm.conn[2,3] = fac[j]*conn[2,3] 
                nm.conn[3,1] = fac[k]*conn[3,1] 
                nm.conn[4,1] = fac[l]*conn[4,1] 
                solve(nm,saveat=nm.hp.saveat, reltol=1e-7, abstol=1e-7)
                println("Status ($(i)/$(length(fac))), ($(j)/$(length(fac))), ($(k)/$(length(fac))), ($(l)/$(length(fac)))")
                plot_analysis(nm,"TC-$(i)-$(j)-$(k)-$(l)")
                plot_syn(nm,1,"time (ms.)","TC-$(i)-$(j)-$(k)-$(l)-Syn")
            end
        end
    end
end


perc_ = [0.85]
# perc_ = [0.8]
vATP_th = [1.1; 1.2; 1.3; 1.4; 1.5; 1.6]
NKA_th = [1.1; 1.2; 1.3; 1.4; 1.5; 1.6]
nm.hp.bandpass = [0.1,40.0]

nm.hp.tstart = 0.5*60*1e3 # ED start time
nm.hp.tend =  10*60*1e3 # ED end time
nm.hp.tfinal = 15*60*1e3 # Simulation end time
nm.hp.saveat = 4
# nm(similar(nm.hp.X0),nm.hp.X0,0,0)

for i in 1:length(perc_)
    for j in 1:length(vATP_th)
        for k in 1:length(NKA_th)
            if vATP_th[j]!=NKA_th[k]
                continue
            else
                nm.hp.perc = perc_[i]
                nm.hp.O2e_th_vATP = vATP_th[j]
                nm.hp.O2e_th_NKA = NKA_th[k]
                solve(nm,saveat=nm.hp.saveat, reltol=1e-7, abstol=1e-7)
                println("Perc ($(i)/$(length(perc_))), vATP($(j)/$(length(vATP_th))), NKA_th($(k)/$(length(NKA_th)))")
                println("Membrane potential: $(nm(:Cortex,"VE")[end])")
                plot_analysis(nm,"Inf-Perc$(Int(perc_[i]*100))-vATP$(vATP_th[j])-NKA$(NKA_th[k])")
                plot_syn(nm,1,"time (ms.)","Inf-Perc$(Int(perc_[i]*100))-vATP$(vATP_th[j])-NKA$(NKA_th[k])-Syn")
            end
        end
    end
end

# perc_ = [0.9; 0.8; 0.7; 0.6; 0.5]
# perc_ = [0.8]
# vATP_th = [1.1; 1.2; 1.3; 1.4; 1.5; 1.6]
# NKA_th = [1.1; 1.2; 1.3; 1.4; 1.5; 1.6]
# nm.hp.bandpass = [0.1,40.0]

# nm(similar(nm.hp.X0),nm.hp.X0,0,0)

# for i in 1:length(perc_)
#     for j in 1:length(vATP_th)
#         for k in 1:length(NKA_th)
#             if vATP_th[j]<NKA_th[k]
#                 continue
#             else
#                 nm.hp.perc = perc_[i]
#                 nm.hp.O2e_th_vATP = vATP_th[j]
#                 nm.hp.O2e_th_NKA = NKA_th[k]
#                 solve(nm,BS3(),reltol=1e-7,abstol=1e-7,saveat=nm.hp.saveat)
#                 println("Perc ($(i)/$(length(perc_))), vATP($(j)/$(length(vATP_th))), NKA_th($(k)/$(length(NKA_th)))")
#                 println("Membrane potential: $(nm(:Cortex,"VE")[end])")
#                 plot_analysis(nm,"Perc$(Int(perc_[i]*100))-vATP$(vATP_th[j])-NKA$(NKA_th[k])")
#             end
#         end
#     end
# end


# perc_ = [0.9; 0.8; 0.7; 0.6; 0.5]
# fac_e = [1.0; 1.5; 2/0; 0.1; 0.5]
# fac_i = [1.0; 1.5; 2.0; 0.1; 0.5]
# nm.hp.bandpass = [0.1,40.0]
# 
# # nm(similar(nm.hp.X0),nm.hp.X0,0,0)
# 
# 
# syn_act_1 = nm.areas[2].pop1.syn_act
# syn_act_2 = nm.areas[2].pop2.syn_act
# syn_deact_1 = nm.areas[2].pop1.syn_deact
# syn_deact_2 = nm.areas[2].pop2.syn_deact
# 
# for i in 1:length(perc_)
#     for j in 1:length(fac_e)
#         for k in 1:length(fac_i)
#                 nm.hp.perc = perc_[i]
#                 nm.areas[2].pop1.syn_act = fac_e[j]*syn_act_1
#                 nm.areas[2].pop2.syn_act = fac_i[k]*syn_act_2
#                 nm.areas[2].pop1.syn_deact = fac_e[j]*syn_deact_1
#                 nm.areas[2].pop2.syn_deact = fac_i[k]*syn_deact_2
#                 solve(nm,BS3(),reltol=1e-7,abstol=1e-7,saveat=nm.hp.saveat)
#                 println("Perc ($(i)/$(length(perc_))), FacE($(j)/$(length(fac_e))), FacI($(k)/$(length(fac_i)))")
#                 plot_analysis(nm,"Perc$(Int(perc_[i]*100))-FacE$(fac_e[j])-FacI$(fac_i[k])")
#         end
#     end
# end


# solve(nm,saveat=hp.saveat,reltol=1e-9,abstol=1e-9) # Solve system using solve(nm)
# solve(nm,CVODE_BDF(),saveat=hp.saveat,reltol=1e-7,abstol=1e-7) # Solve system using solve(nm)
# plot_syn(nm,1,"time (ms.)","ModerateED") # Plot synaptic currents, change second argument to 60*1000 to plot in min.
# plot(nm,"Baseline") # Plots EEG and ion dynamics for both regions and saves as two files.
