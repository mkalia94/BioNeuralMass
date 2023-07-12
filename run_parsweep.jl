using Pkg
dir = pwd()
Pkg.activate("$(dir)/.")
using BioNeuralMass, Sundials, DifferentialEquations, Plots, Waveforms, DSP, SignalAnalysis, Peaks

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
nm.conn = [ 0.0 2  0.5 0.0 ; 0.3 0.0 0.1 0.0; 0.5 0.0  0.3 5; 0.2 0 0.5 2.5] #Thalamus first

# Thalamus baseline conditions 
thalamus_.pop1.syn_th = 0.2 
thalamus_.pop1.syn_act = 1.25
thalamus_.pop1.syn_deact = 0.3
thalamus_.pop2.syn_th = 0.5
thalamus_.pop2.syn_act = 0.5
thalamus_.pop2.syn_deact = 0.003


# Hyperparameters
hp.excite = [20, 1000,20*60*1e3, 1e8] # [Current strength, start time, end time, large number]
hp.perc = 0.7# Available energy ∈ [0,1], 0 -> complete ED, 1 -> No ED
hp.tstart = 1*60*1e3 # ED start time
hp.tend =  3*60*1e3 # ED end time
hp.tfinal = 5*60*1e3 # Simulation end time
hp.beta1 = 5 # ED onset rate (higher is faster)
hp.beta2 = 5 # ED offset rate (higher is faster)
hp.saveat = 4 # save every at every x milliseconds
hp.O2e_th_vATP = 1.5 # vATPase oxygen threshold
hp.O2e_th_NKA = 1.5 # NKA oxygen threshold

nm.areas = [thalamus_,cortex_] # Order of areas matters! First area is the one stimulated. Take care of nm.conn!
nm.hp = hp


###################### Effect of stimulation input on EEG  #########################
# Stim = range(41,60,length=20)
# freqs = zeros(length(Stim))
# amps = zeros(length(Stim))
# nm.hp.excite = [20, 1000, 1e8, 1e8]
# nm.hp.perc = 1
# nm.hp.saveat = 0.5
# nm.hp.tfinal = 1e4

# for i in 1:length(Stim)
#    nm.hp.excite[1] = Stim[i]
#    solve(nm,saveat=nm.hp.saveat, reltol=1e-7, abstol=1e-7)
#    EEG = nm("EEGraw")[:,1]
#    peaks, _ = findmaxima(EEG)
#    # val, ind = findmax(diff(peaks))
#    val = diff(peaks)[end]
#    freqs[i] = 1000/(val*nm.hp.saveat)
#    amps[i] = maximum(EEG) - minimum(EEG)
#    println("$(i) done, frequency: $(freqs[i])")
# end

###################################################################################


######################## Effect of oxygen thresholds - two parameter sweep ##########

# perc_ = [0.8]
# vATP_th = [1.1; 1.2; 1.3; 1.4; 1.5; 1.6]
# NKA_th = [1.1; 1.2; 1.3; 1.4; 1.5; 1.6]
# nm.hp.bandpass = [0.1,40.0]

# nm.hp.tstart = 0.5*60*1e3 ED start time
# nm.hp.tend =  10*60*1e3 ED end time
# nm.hp.tfinal = 15*60*1e3 Simulation end time
# nm.hp.saveat = 4
# nm(similar(nm.hp.X0),nm.hp.X0,0,0)

# for i in 1:length(perc_)
#     for j in 1:length(vATP_th)
#         for k in 1:length(NKA_th)
#             if vATP_th[j] < NKA_th[k]
#                 continue
#             else
#                 nm.hp.perc = perc_[i]
#                 nm.hp.O2e_th_vATP = vATP_th[j]
#                 nm.hp.O2e_th_NKA = NKA_th[k]
#                 solve(nm,saveat=nm.hp.saveat)
#                 println("Perc ($(i)/$(length(perc_))), vATP($(j)/$(length(vATP_th))), NKA_th($(k)/$(length(NKA_th)))")
#                 println("Membrane potential: $(nm(:Cortex,"VE")[end])")
#                 plot_analysis(nm,"Inf-Perc$(Int(perc_[i]*100))-vATP$(vATP_th[j])-NKA$(NKA_th[k])")
#                 plot_syn(nm,1,"time (ms.)","Inf-Perc$(Int(perc_[i]*100))-vATP$(vATP_th[j])-NKA$(NKA_th[k])-Syn")
#             end
#         end
#     end
# end


# #########################################################################################


# #################### Three examples from parameter sweep ##############################
perc_ = [0.8]
nm.hp.bandpass = [0.1,40.0]

nm.hp.tstart = 0.5*60*1e3 # ED start time
nm.hp.tend =  10*60*1e3 # ED end time
nm.hp.tfinal = 15*60*1e3 # Simulation end time
nm.hp.saveat = 1

th = [1.5 1.1; 1.3 1.3; 1.6 1.6] 

for i in 1:length(perc_)
    for j in [1,2,3]
                nm.hp.perc = perc_[i]
                nm.hp.O2e_th_vATP = th[j,1]
                nm.hp.O2e_th_NKA = th[j,2]
                solve(nm,saveat=nm.hp.saveat,reltol=1e-6,abstol=1e-6)
                responsetype = Bandpass(nm.hp.bandpass[1],nm.hp.bandpass[2],fs=1000/(nm.hp.saveat))
                designmethod = Butterworth(4)
                EEG = filt(digitalfilter(responsetype,designmethod),nm("EEGraw"))
                EEG = reshape(EEG,length(EEG))
                y = tfd(EEG,Spectrogram(nfft=Int(3000/nm.hp.saveat),noverlap=Int(1000/nm.hp.saveat),window=DSP.hamming); fs=Int(1000/nm.hp.saveat))
                ind = findmax(y.power,dims=1)[2]
                freqs = [y.freq[ind[1,i][1]] for i in 1:length(y.time)]
                plot(y.time/60/1e3,freqs,label="",lw=2,c=:black,ylims=(0,13))
                plot!(size=(300,200))
                savefig("freq_$(j).pdf")
                plot(nm.sol.t/60/1e3,EEG,label="",lw=2,c=:black)
                plot!(size=(400,200))
                savefig("EEG_$(j).pdf")
               VCE = nm(:Cortex,"VE")
               VCI = nm(:Cortex,"VI")
               plot(nm.sol.t/60/1e3,VCE,c=:blue,lw=2,label="")
               plot!(nm.sol.t/60/1e3,VCI,c=:red,lw=2,label="")
               plot!(size=(300,200))
               savefig("MemPot_$(j).pdf")
               println("$(j) done")
   end
end

# ######################################################################################


# ####################### Effect of thalamocortical input - two parameter sweep #######


# fac2 = range(0,2,length=50)
# fac1 = range(0.2,1.2,length=25)
# freqs = zeros(25,50)
# amps = zeros(25,50)
# conn = deepcopy(nm.conn)

# for i in 1:length(fac1)
#     for j in 1:length(fac2)
#         nm.hp.perc = 1
#         nm.hp.tfinal = 5e3
#         nm.hp.saveat = 0.1
#         nm.conn[3,1] = fac1[i]*conn[3,1] 
#         nm.conn[4,1] = fac2[j]*conn[4,1] 
#         solve(nm,saveat=nm.hp.saveat, reltol=1e-7, abstol=1e-7)
#         println("Status ($(i)/$(length(fac))), ($(j)/$(length(fac)))")
#         EEG = nm("EEGraw")[:,1]
#         pp = periodogram(EEG,fs = 1000/nm.hp.saveat)
#         peaks,vals = findmaxima(pp.power)
#         _, ind = findmax(vals)
#         if pp.freq[peaks[ind]] > 40
#             continue
#         end
#         freqs[i,j] = pp.freq[peaks[ind]]
#         amps[i,j] = maximum(EEG[end-1000:end])-minimum(EEG[end-1000:end])
#     end
# end
# xs = [string("$(i)"[1:4]) for i in 0.500001 .* fac1]
# ys = ["0.0"; [string("$(i)"[1:4]) for i in (0.2001 .* fac2)[2:end]]]
# heatmap(xs,ys,freq',c=:grayC,tickfontsize=12,right_margin=3Plots.mm,reuse=false)

# ######################################################################################



# ####################### Left MMOs - parameter sweep ########################## #######

# npeaks = zeros(100)
# periods = zeros(100)
# conn = deepcopy(nm.conn)
# conn41 = range(0.11,0.13,length=100)
# EEG = nothing


# for i in [90,91]
#     global EEG
#     nm.hp.perc = 1
#     nm.hp.tfinal = 1e4
#     nm.hp.saveat = 0.1
#     nm.conn[3,1] = 0.225 
#     nm.conn[4,1] = conn41[i]
#     try
#         solve(nm,saveat=nm.hp.saveat, reltol=1e-7, abstol=1e-7)
#     catch e
#         continue
#     end
#     EEG = nm("EEGraw")[:,1]
#     peaks, _ = findmaxima(EEG[80000:100000])
#      val, ind = findmax(diff(peaks))
#     xx = findall(x->abs(x-val)/val<0.1,diff(peaks))
#     try 
#         npeaks[i] = diff(xx)[end]
#     catch e
#          continue
#     end
#     if ind - diff(xx)[end]<1
#         periods[i] = sum(diff(peaks)[ind+1:ind + diff(xx)[end]])
#     else
#         periods[i] = sum(diff(peaks)[ind-diff(xx)[end]+1:ind])
#     end
#     println("Status ($(i)/$(length(conn41))), nPeaks: $(npeaks[i]), Period: $(periods[i])")
# end#

# ######################################################################################



# ####################### Left MMOs - parameter sweep ########################## #######

# npeaks2 = zeros(100)
# periods2 = zeros(100)
# conn = deepcopy(nm.conn)
# conn41 = range(0,0.04,length=100)
#  
# for i in 1:length(conn41)
#     nm.hp.perc = 1
#     nm.hp.tfinal = 1e5
#     nm.hp.saveat = 0.1
#     nm.conn[3,1] = 0.475 
#     nm.conn[4,1] = conn41[i]
#     try
#         solve(nm,saveat=nm.hp.saveat, reltol=1e-7, abstol=1e-7)
#     catch e
#         continue
#     end
#     EEG = nm("EEGraw")[:,1]
#     l = length(EEG)
#     peaks, _ = findmaxima(EEG[l-20000:l])
#     val, ind = findmax(diff(peaks))
#     xx = findall(x->abs(x-val)/val<0.1,diff(peaks))
#     try 
#         npeaks2[i] = diff(xx)[end]
#         if ind - diff(xx)[end]<1
#            periods2[i] = sum(diff(peaks)[ind+1:ind + diff(xx)[end]])
#         else
#            periods2[i] = sum(diff(peaks)[ind-diff(xx)[end]+1:ind])
#         end
#     catch e
#         npeaks2[i] = xx[end]
#         if ind - xx[end]<1
#            periods2[i] = sum(diff(peaks)[ind+1:ind + xx[end]])
#         else
#            periods2[i] = sum(diff(peaks)[ind-xx[end]+1:ind])
#         end
#     end
#          println("Status ($(i)/$(length(conn41))), nPeaks: $(npeaks2[i]), Period: $(periods2[i])")
# end
 





