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
nm.conn = [5.0  3.0; 4.0 0.5]
hp.X0 = [pop1.NaCi0*pop1.Wi0,pop1.KCi0*pop1.Wi0,pop1.ClCi0*pop1.Wi0,pop1.Wi0]
hp.X0 = [hp.X0; hp.X0; 0.008;0.008]

hp.excite = [20.0, 1000.0,1*60*1000.0 + 2*60*1000, 1e10]
hp.tstart = 1000
hp.tend = 1*60*1000 + 2*1000*60
hp.tfinal = 4*60*1000
hp.beta1 = 0.004
hp.beta2 = 0.004
hp.saveat = 5

perc = [1]

for p in perc
    hp.perc = p
    try
        solve(nm,hp,[cortex_,],reltol=1e-9,abstol=1e-9)
    catch e
    end
    if !ismissing(nm.sol)
        plot(hp,nm,[cortex_,],"Perc$(p*100)")
    end
end

