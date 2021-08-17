using Pkg
dir = pwd()
Pkg.activate("$(dir)/.")
using BioNeuralMass, Sundials, DifferentialEquations, Plots, Waveforms

# Set up cortical neural area
hp = HyperParam(synapseoff = false, ratio = 0.8)
pop1 = NeuralPopSoma{Cortex, Excitatory}()
pop2 = NeuralPopSoma{Cortex, Inhibitory}()
cortex_ = NeuralArea(hp,pop1,pop2);
Par(hp,cortex_)

# Set up neural mass
nm = BioNM()
nm.conn = [0.7  5;0.5  2.5]

hp.excite = [15.0, 1000.0,1000.0 + 2*60*1000, 1e10]
hp.tstart = 1000
hp.tend = 1000+2*60*1000
hp.tfinal = 3*60*1000
hp.beta1 = 4
hp.beta2 = 4
hp.saveat = 1

hp.perc = 0.7

solve(nm,hp,[cortex_,],saveat=hp.saveat,reltol=1e-8,abstol=1e-8)
syn_curr = nm([cortex_,],hp,"syn_curr")
plot(nm.sol.t/60/1000,syn_curr,layout = (2,1))

# for p in perc
#     hp.perc = p
#     try
#         solve(nm,hp,[cortex_,],reltol=1e-9,abstol=1e-9)
#     catch e
#     end
#     if !ismissing(nm.sol)
#         plot(hp,nm,[cortex_,],"Perc$(p*100)")
#     end
# end

