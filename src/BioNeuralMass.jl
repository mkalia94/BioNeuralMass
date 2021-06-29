module BioNeuralMass

using DifferentialEquations, Plots, QuadGK, LaTeXStrings, DSP, Waveforms
include("build.jl")
include("model.jl")
include("params.jl")
include("utils.jl")

export HyperParam, InitVal, BioNM
export Par
export ampa, gaba, FiringRate
export solve, plot

end # module
