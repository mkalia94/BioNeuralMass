using Pkg
dir = pwd()
Pkg.activate("$(dir)/.")
using BioNeuralMass, Sundials, DifferentialEquations, Plots, Waveforms

hp_ = HyperParam(synapseoff=false)
in_ = InitVal()
nm = BioNM(hp = hp_, init = in_)
Par(nm)

nm.hp.excite = [20.0, 1000.0,1000.0 + 2*60*1000, 1e10]
nm.hp.perc = 1
nm.hp.tstart = 1000
nm.hp.tend = 1000 + 2*60*1000
nm.hp.tfinal = 5*60*1000
nm.hp.beta1 = 0.004
nm.hp.beta2 = 0.004

nm.g_EE = nm.g_EE * 1e5
nm.g_EI = nm.g_EI * 1e5
nm.g_II = nm.g_II * 1e5
nm.g_IE = nm.g_IE * 1e5

nm.g_IE = 3
nm.g_EE = 5
nm.g_EI = 4

nm.b_ampa = 0.05
nm.b_gaba = 0.05

nm.ampa_th = 0.1
nm.gaba_th = 0.5

nm.hp.saveat = 1

perc_ = [1, 0.8, 0.6, 0.4, 0.2]

for p in perc_
    nm.hp.perc = p
    try 
        solve(nm,reltol=1e-10,abstol=1e-10,saveat=nm.hp.saveat)
    catch e
    end  
    plot(nm)
    savefig("Perc$(p*100).pdf")
end
