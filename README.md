# BioNeuralMass

Please see `run.jl` in order to run the simulations shown in the Overleaf document. Here's an overview of the code:

1. `NeuralPopSoma{Area,Type}` is a method that creates a point neuron model in region `Area' (can be either `Cortex` or `Thalamus` at this point) of type `Type` (`Inhibitory/Excitatory`)
2. Two `NeuralPopSoma` can be grouped into one `NeuralArea`, for which the function `Par` computes all impermeants and leak conductances
3. The method `BioNM` constructs the neural mass that couples multiple areas, via `BioNM.conn`. Note that `BioNM.areas` is a vector of type `NeuralArea`. The first element of this vector receives stimulation
4. Change the hyperparameters using the method `HyperParam`
5. `solve(nm::BioNM)`, `plot(nm::BioNM)`, `plot_syn(nm::BioNM,1)` solve the system, plot ion dynamics and plot synaptic currents respectively.
