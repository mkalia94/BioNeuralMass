# Define hyper parameters
Base.@kwdef mutable struct HyperParam
    tstart ::Union{Float64,Missing} = 20.0 
    tend   ::Union{Float64,Missing} = 30.0 
    perc   ::Union{Float64,Missing} = 0.5 
    beta1  ::Union{Float64,Missing} = 4  
    beta2  ::Union{Float64,Missing} = 4
    excite :: Union{Array{Float64},Missing} = missing 
    tfinal :: Union{Float64,Missing} = 50.0 
    ratio  ::Union{Float64,Missing} = missing
    saveat :: Union{Float64,Missing} = 0.01
    bandpass :: Union{Array{Float64},Missing} = [5.0; 40.0]
    synapseoff :: Bool
end             

# Define initial conditions
Base.@kwdef mutable struct InitVal
    NaCiE :: Union{Float64,Missing} = missing
    NaCiI :: Union{Float64,Missing} = missing
    KCiE :: Union{Float64,Missing} = missing
    KCiI :: Union{Float64,Missing} = missing
    ClCiE :: Union{Float64,Missing} = missing
    ClCiI :: Union{Float64,Missing} = missing
    rAMPA  :: Union{Float64,Missing} = missing
    rGABA  :: Union{Float64,Missing} = missing
end

# Define all model parameters and main method
Base.@kwdef mutable struct BioNM
        
    rcell         ::Union{Float64,Missing}    = 7.815926417967719            # [mu m], radius of spherical cell
    Acell         ::Union{Float64,Missing}    = missing                       # [mu m^2], Membrane surface of neuron
    Wcell         ::Union{Float64,Missing}    = missing                       # [mu m^3], Volume of spherical cell
    Ari0  ::Union{Float64,Missing}    = missing                       # [1000 mu m^2]
    Wi0  ::Union{Float64,Missing}    = missing                       # [1000 mu m^3]
    We0  ::Union{Float64,Missing}    = missing                       # [1000 mu m^3]
    Wtot ::Union{Float64,Missing} = missing
    C             ::Union{Float64,Missing}    = 20                           # [pF], membrane capacitance
    F             ::Union{Float64,Missing}    = 96485.333                    # [C/mol], Faraday's constant
    R             ::Union{Float64,Missing}    = Float64(8.3144598*1000)               # [(CmV)/(mol K)], Universal gas constant
    T             ::Union{Float64,Missing}    = 310                          # [K], Absolute temperature
    xi            ::Union{Float64,Missing}    = missing                       # [1/mV]

    # Part 1B. ATP-related parameters and ATP-pump
    rho_p_E       ::Union{Float64,Missing}    = 60                                 # [pA], Pump current scaling (including strenght of pump) / NaK pump rate
    rho_p_I       ::Union{Float64,Missing}    = 60                                 # [pA], Pump current scaling (including strenght of pump) / NaK pump rate
    nka_na        ::Union{Float64,Missing}    = 13                            
    nka_k        ::Union{Float64,Missing}    =0.2                             

    
    # Part 1C. Transmembrane receptors
    a_ampa     ::Union{Float64,Missing}   = 2          # [1/ms], forward rate constant                                   
    b_ampa     ::Union{Float64,Missing}   = 0.5           # [1/ms], backward rate constant
    g_EE   ::Union{Float64,Missing}   = 2f-5 #Float64(0.1*1.7)         # [nS], maximal conductance, o.1 - 1 nS (pS is voor 1 single channel
    g_EI   ::Union{Float64,Missing}   = 1f-5                                                                            
    a_gaba     ::Union{Float64,Missing}   = 14          # [1/ms]
    b_gaba     ::Union{Float64,Missing}   = 0.03         # 1/ms  --> slower than nm.b_ampa, because gaba should be slower!
    g_IE   ::Union{Float64,Missing}   = 1f-5  
    g_II   ::Union{Float64,Missing}   = 5f-6          # [nS] 0.25 - 1.2 nS 
    ampa_th ::Union{Float64,Missing} = 0.5
    gaba_th ::Union{Float64,Missing} = 0.9 

    # Part 1D. Connection parameters 1/(nm.F)
    N_ee_b     ::Union{Float64,Missing}   = 3000    # [-], number of synaptic contacts from excitatory to exctitatory neurons
    N_ei_b     ::Union{Float64,Missing}   = 3000    # [-], number of synaptic contacts from excitatory to inhibitory neurons
    N_ie_b     ::Union{Float64,Missing}   = 500     # [-], number of synaptic contacts from inhibitory to exctitatory neurons
    N_ii_b     ::Union{Float64,Missing}   = 500     # [-], number of synaptic contacts from inhibitory to inhibitory neurons

    # Part 1E. KCL cotransporter
    U_kcl      ::Union{Float64,Missing} = 1f-6 # [fmol/ms/mV] the cotransporter strength 
    PWi        ::Union{Float64,Missing} = 2*1f-14   

    # Part 1F. Ion related parameters
    NaCe0     ::Union{Float64,Missing} = 152   # [mMol/l = mM], Initial concentration of Sodium (o = extracellular, e=excitatory) 
    NaCiE0   ::Union{Float64,Missing} = 13    # [mMol/l = mM], ICS sodium concentration of excitatory population
    NaCiI0   ::Union{Float64,Missing} = 13    # [mMol/l = mM], ICS sodium concentration of inhibitory population
                                                                                                      
    KCe0      ::Union{Float64,Missing} = 3    # [mMol/l = mM], ECS potassium concentration
    KCiE0    ::Union{Float64,Missing} = 145   # [mMol/l = mM], Initial ICS potassium concentration of excitatory population
    KCiI0    ::Union{Float64,Missing} = 145   # [mMol/l = mM], Initial ICS potassium concentration of inhibitory population
    ClCe0     ::Union{Float64,Missing} = 135                                                                        
    ClCiE0   ::Union{Float64,Missing} = 8.0   # [mMol/l = mM], ICS chloride concentration of excitatory population
    ClCiI0   ::Union{Float64,Missing} = 8.0   # [mMol/l = mM], ICS chloride concentration of inhibitory population

    Vi0 ::Union{Float64,Missing} = -65.5 
                                             
    CNa       ::Union{Float64,Missing} = missing  
    CK        ::Union{Float64,Missing} = missing  
    CCl       ::Union{Float64,Missing} = missing  
    E_na_E_0   ::Union{Float64,Missing} = missing  # [mV], Initial Nernst potential of sodium                             
    E_na_I_0   ::Union{Float64,Missing} = missing  # [mV], Initial Nernst potential of sodium of inhibitory population
    
    E_k_E_0    ::Union{Float64,Missing} = missing  # [mV], Initial Nernst potential of potassium
    E_k_I_0    ::Union{Float64,Missing} = missing  # [mV], Initial Nernst potential of potassium of inhibitory population

    # Part 1G. Leak permeabilities (with voltage-gated currents)
    PnaL       ::Union{Float64,Missing}  = 0.001131718785164   # [1000 mum^3/ms], Leak Na+ permeability (for -65 mV) 
    PkL        ::Union{Float64,Missing}  = 0.013554843976053   # [1000 mum^3/ms], Leak K+ permeability (for -65mV)
    PclL       ::Union{Float64,Missing}  = 0.001930766455246   # [1000 mum^3/ms], Leak Cl- permeability (for -65mV)

    # Part 1H. Permeabilities
    PnaT       ::Union{Float64,Missing} = 80*1f-5                               # [1000 mum^3/ms] Maximal transient Na+ permeability         
    PkD        ::Union{Float64,Missing} = 40*1f-5                               # [1000 mum^3/ms] Maximal delayed rectifier K+ permeability 
    PclG       ::Union{Float64,Missing} = 1.95*1f-5                             # [1000 mum^3/ms] Maximal voltage-gated Cl- permeability     

    # Part 1I. Firing Rate related parameters
    sigma      ::Union{Float64,Missing} = 1 

    # Part 1Ia. Coefficients of polynomial Ithreshold
    p00        ::Union{Float64,Missing}       = -213.8 
    p10        ::Union{Float64,Missing}      = -8.404 
    p01        ::Union{Float64,Missing}      = -3.007
    p20        ::Union{Float64,Missing}      = -0.1111 
    p11        ::Union{Float64,Missing}      = -0.07862
    p02        ::Union{Float64,Missing}      = -0.001213
    p30        ::Union{Float64,Missing}      = -0.0005002 
    p21        ::Union{Float64,Missing}      = -0.0005221 
    p12        ::Union{Float64,Missing}      = -1.483f-06 
    p03        ::Union{Float64,Missing}      = 3.779f-05 

    # Part 1Ib. Coefficients of polynomial 'kappa'
    q00       ::Union{Float64,Missing} = 0.287 
    q10       ::Union{Float64,Missing} = 0.01094
    q01       ::Union{Float64,Missing} = 0.001027
    q20       ::Union{Float64,Missing} = 0.0001692
    q11       ::Union{Float64,Missing} = -0.00000842
    q02       ::Union{Float64,Missing} = -0.00004839
    q30       ::Union{Float64,Missing} = 0.0000008167
    q21       ::Union{Float64,Missing} = -0.0000003156
    q12       ::Union{Float64,Missing} = -0.0000004421
    q03       ::Union{Float64,Missing} = 0.0000002764

    # Part 1Ic. Coefficients of polynomial Ithreshold2
    r00      ::Union{Float64,Missing}   = -185.9
    r10      ::Union{Float64,Missing}   = -4.877
    r01      ::Union{Float64,Missing}   = -3.674
    r20      ::Union{Float64,Missing}   = -0.002785
    r11      ::Union{Float64,Missing}   = -0.1044
    r02      ::Union{Float64,Missing}   = 0.01561
    r30      ::Union{Float64,Missing}   = -1.17f-05
    r21      ::Union{Float64,Missing}   = 1.897f-05
    r12      ::Union{Float64,Missing}   = 0.0006927
    r03      ::Union{Float64,Missing}   = 0.0001089

    NAe     :: Union{Float64,Missing} = missing 
    NAi     :: Union{Float64,Missing} = missing

    hp :: Union{HyperParam,Missing} = missing
    init :: Union{InitVal,Missing} = missing
    
    X0 :: Union{Array{Float64},Missing} = missing
    sol :: Union{ODESolution,Missing} = missing

end

