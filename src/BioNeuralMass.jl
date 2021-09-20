module BioNeuralMass

using DifferentialEquations, Plots, QuadGK, LaTeXStrings, DSP, Waveforms, Quadrature
include("pop.jl")
include("neuralarea.jl")
include("params.jl")
include("synapse.jl")
include("utils.jl")


export Area, Behaviour, Excitatory, Inhibitory, Thalamus, Cortex
export HyperParam, NeuralPopSoma
export NeuralArea, get_ECS
export Par 
export BioNM
export solve, plot,plot_syn


end # module
