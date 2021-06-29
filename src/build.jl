# Define hyper parameters
Base.@kwdef mutable struct HyperParam
    tstart ::Union{Float32,Missing} = 20.0f0 
    tend   ::Union{Float32,Missing} = 30.0f0 
    perc   ::Union{Float32,Missing} = 0.5f0 
    beta1  ::Union{Float32,Missing} = 4f0  
    beta2  ::Union{Float32,Missing} = 4f0
    excite :: Union{Array{Float32},Missing} = missing 
    tfinal :: Union{Float32,Missing} = 50.0f0 
    ratio  ::Union{Float32,Missing} = missing
    saveat :: Union{Float32,Missing} = 0.01f0
    bandpass :: Union{Array{Float32},Missing} = [5.0f0; 40.0f0]
end             

# Define initial conditions
Base.@kwdef mutable struct InitVal
    NaCiE :: Union{Float32,Missing} = missing
    NaCiI :: Union{Float32,Missing} = missing
    KCiE :: Union{Float32,Missing} = missing
    KCiI :: Union{Float32,Missing} = missing
    ClCiE :: Union{Float32,Missing} = missing
    ClCiI :: Union{Float32,Missing} = missing
    rAMPA  :: Union{Float32,Missing} = missing
    rGABA  :: Union{Float32,Missing} = missing
end

# Define all model parameters and main method
Base.@kwdef mutable struct BioNM
        
    rcell         ::Union{Float32,Missing}    = 7.815926417967719f0            # [mu m], radius of spherical cell
    Acell         ::Union{Float32,Missing}    = missing                       # [mu m^2], Membrane surface of neuron
    Wcell         ::Union{Float32,Missing}    = missing                       # [mu m^3], Volume of spherical cell
    Ari0  ::Union{Float32,Missing}    = missing                       # [1000 mu m^2]
    Wi0  ::Union{Float32,Missing}    = missing                       # [1000 mu m^3]
    We0  ::Union{Float32,Missing}    = missing                       # [1000 mu m^3]
    Wtot ::Union{Float32,Missing} = missing
    C             ::Union{Float32,Missing}    = 20f0                           # [pF], membrane capacitance
    F             ::Union{Float32,Missing}    = 96485.333f0                    # [C/mol], Faraday's constant
    R             ::Union{Float32,Missing}    = Float32(8.3144598*1000)               # [(CmV)/(mol K)], Universal gas constant
    T             ::Union{Float32,Missing}    = 310f0                          # [K], Absolute temperature
    xi            ::Union{Float32,Missing}    = missing                       # [1/mV]

    # Part 1B. ATP-related parameters and ATP-pump
    rho_p_E       ::Union{Float32,Missing}    = 60                                 # [pA], Pump current scaling (including strenght of pump) / NaK pump rate
    rho_p_I       ::Union{Float32,Missing}    = 60                                 # [pA], Pump current scaling (including strenght of pump) / NaK pump rate
    nka_na        ::Union{Float32,Missing}    = 13                            
    nka_k        ::Union{Float32,Missing}    =0.2                             

    
    # Part 1C. Transmembrane receptors
    a_ampa     ::Union{Float32,Missing}   = 2f0          # [1/ms], forward rate constant                                   
    b_ampa     ::Union{Float32,Missing}   = 0.5f0           # [1/ms], backward rate constant
    g_ampa_E   ::Union{Float32,Missing}   = Float32(0.1*1.7)         # [nS], maximal conductance, o.1 - 1 nS (pS is voor 1 single channel
    g_ampa_I   ::Union{Float32,Missing}   = 0.1f0                                                                            
    a_gaba     ::Union{Float32,Missing}   = 14f0          # [1/ms]
    b_gaba     ::Union{Float32,Missing}   = 0.03f0         # 1/ms  --> slower than nm.b_ampa, because gaba should be slower!
    g_gaba_E   ::Union{Float32,Missing}   = 0.1f0  
    g_gaba_I   ::Union{Float32,Missing}   = 0.05f0          # [nS] 0.25 - 1.2 nS  

    # Part 1D. Connection parameters
    N_ee_b     ::Union{Float32,Missing}   = 3000f0    # [-], number of synaptic contacts from excitatory to exctitatory neurons
    N_ei_b     ::Union{Float32,Missing}   = 3000f0    # [-], number of synaptic contacts from excitatory to inhibitory neurons
    N_ie_b     ::Union{Float32,Missing}   = 500f0     # [-], number of synaptic contacts from inhibitory to exctitatory neurons
    N_ii_b     ::Union{Float32,Missing}   = 500f0     # [-], number of synaptic contacts from inhibitory to inhibitory neurons

    # Part 1E. KCL cotransporter
    U_kcl      ::Union{Float32,Missing} = 1f-6 # [fmol/ms/mV] the cotransporter strength 
    PWi        ::Union{Float32,Missing} = 2*1f-14   

    # Part 1F. Ion related parameters
    NaCe0     ::Union{Float32,Missing} = 152f0   # [mMol/l = mM], Initial concentration of Sodium (o = extracellular, e=excitatory) 
    NaCiE0   ::Union{Float32,Missing} = 13f0    # [mMol/l = mM], ICS sodium concentration of excitatory population
    NaCiI0   ::Union{Float32,Missing} = 13f0    # [mMol/l = mM], ICS sodium concentration of inhibitory population
                                                                                                      
    KCe0      ::Union{Float32,Missing} = 3f0    # [mMol/l = mM], ECS potassium concentration
    KCiE0    ::Union{Float32,Missing} = 145f0   # [mMol/l = mM], Initial ICS potassium concentration of excitatory population
    KCiI0    ::Union{Float32,Missing} = 145f0   # [mMol/l = mM], Initial ICS potassium concentration of inhibitory population
    ClCe0     ::Union{Float32,Missing} = 135f0                                                                        
    ClCiE0   ::Union{Float32,Missing} = 8.0f0   # [mMol/l = mM], ICS chloride concentration of excitatory population
    ClCiI0   ::Union{Float32,Missing} = 8.0f0   # [mMol/l = mM], ICS chloride concentration of inhibitory population

    Vi0 ::Union{Float32,Missing} = -65.5f0 
                                             
    CNa       ::Union{Float32,Missing} = missing  
    CK        ::Union{Float32,Missing} = missing  
    CCl       ::Union{Float32,Missing} = missing  
                                             
    E_na_E_0   ::Union{Float32,Missing} = missing  # [mV], Initial Nernst potential of sodium                             
    E_na_I_0   ::Union{Float32,Missing} = missing  # [mV], Initial Nernst potential of sodium of inhibitory population
    
    E_k_E_0    ::Union{Float32,Missing} = missing  # [mV], Initial Nernst potential of potassium
    E_k_I_0    ::Union{Float32,Missing} = missing  # [mV], Initial Nernst potential of potassium of inhibitory population

    # Part 1G. Leak permeabilities (with voltage-gated currents)
    PnaL       ::Union{Float32,Missing}  = 0.001131718785164f0   # [1000 mum^3/ms], Leak Na+ permeability (for -65 mV) 
    PkL        ::Union{Float32,Missing}  = 0.013554843976053f0   # [1000 mum^3/ms], Leak K+ permeability (for -65mV)
    PclL       ::Union{Float32,Missing}  = 0.001930766455246f0   # [1000 mum^3/ms], Leak Cl- permeability (for -65mV)

    # Part 1H. Permeabilities
    PnaT       ::Union{Float32,Missing} = 80*1f-5                               # [1000 mum^3/ms] Maximal transient Na+ permeability         
    PkD        ::Union{Float32,Missing} = 40*1f-5                               # [1000 mum^3/ms] Maximal delayed rectifier K+ permeability 
    PclG       ::Union{Float32,Missing} = 1.95*1f-5                             # [1000 mum^3/ms] Maximal voltage-gated Cl- permeability     

    # Part 1I. Firing Rate related parameters
    sigma      ::Union{Float32,Missing} = 1f0 

    # Part 1Ia. Coefficients of polynomial Ithreshold
    p00        ::Union{Float32,Missing}       = -213.8f0 
    p10        ::Union{Float32,Missing}      = -8.404f0 
    p01        ::Union{Float32,Missing}      = -3.007f0
    p20        ::Union{Float32,Missing}      = -0.1111f0 
    p11        ::Union{Float32,Missing}      = -0.07862f0
    p02        ::Union{Float32,Missing}      = -0.001213f0
    p30        ::Union{Float32,Missing}      = -0.0005002f0 
    p21        ::Union{Float32,Missing}      = -0.0005221f0 
    p12        ::Union{Float32,Missing}      = -1.483f-06 
    p03        ::Union{Float32,Missing}      = 3.779f-05 

    # Part 1Ib. Coefficients of polynomial 'kappa'
    q00       ::Union{Float32,Missing} = 0.287f0 
    q10       ::Union{Float32,Missing} = 0.01094f0
    q01       ::Union{Float32,Missing} = 0.001027f0
    q20       ::Union{Float32,Missing} = 0.0001692f0
    q11       ::Union{Float32,Missing} = -0.00000842f0
    q02       ::Union{Float32,Missing} = -0.00004839f0
    q30       ::Union{Float32,Missing} = 0.0000008167f0
    q21       ::Union{Float32,Missing} = -0.0000003156f0
    q12       ::Union{Float32,Missing} = -0.0000004421f0
    q03       ::Union{Float32,Missing} = 0.0000002764f0

    # Part 1Ic. Coefficients of polynomial Ithreshold2
    r00      ::Union{Float32,Missing}   = -185.9f0
    r10      ::Union{Float32,Missing}   = -4.877f0
    r01      ::Union{Float32,Missing}   = -3.674f0
    r20      ::Union{Float32,Missing}   = -0.002785f0
    r11      ::Union{Float32,Missing}   = -0.1044f0
    r02      ::Union{Float32,Missing}   = 0.01561f0
    r30      ::Union{Float32,Missing}   = -1.17f-05
    r21      ::Union{Float32,Missing}   = 1.897f-05
    r12      ::Union{Float32,Missing}   = 0.0006927f0
    r03      ::Union{Float32,Missing}   = 0.0001089f0

    NAe     :: Union{Float32,Missing} = missing 
    NAi     :: Union{Float32,Missing} = missing

    hp :: Union{HyperParam,Missing} = missing
    init :: Union{InitVal,Missing} = missing
    
    X0 :: Union{Array{Float32},Missing} = missing
    sol :: Union{ODESolution,Missing} = missing

end

