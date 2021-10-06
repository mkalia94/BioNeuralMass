abstract type Area end
abstract type Behaviour end
abstract type Excitatory <: Behaviour end
abstract type Inhibitory <: Behaviour end
abstract type Thalamus <: Area end
abstract type Cortex <: Area end

Base.@kwdef mutable struct NeuralPopSoma{A<:Area, B<:Behaviour} 
    C             ::Union{Float64,Missing}    = 20                           # [pF], membrane capacitance
    F             ::Union{Float64,Missing}    = 96485.333                    # [C/mol], Faraday's constant
    R             ::Union{Float64,Missing}    = Float64(8.3144598*1000)               # [(CmV)/(mol K)], Universal gas constant
    T             ::Union{Float64,Missing}    = 310                          # [K], Absolute temperature
    rcell         ::Union{Float64,Missing}    = 7.815926417967719            # [mu m], radius of spherical cell
    Ari0  ::Union{Float64,Missing}    =  4*pi*(rcell)^2/1000;             # [mu m^2], Membrane surface of neuroning                       # [1000 mu m^2]
    Wi0  ::Union{Float64,Missing}    = (4/3)*pi*(rcell)^3/1000;         # [mu m^3], Volume of spherical cell   ng                       # [1000 mu m^3]
    NAi :: Union{Float64,Missing} = missing
    # Part 1B. ATP-related parameters and ATP-pump
    PumpStrength       ::Union{Float64,Missing}    = 60                                 # [pA], Pump current scaling (including strenght of pump) / NaK pump rate
    nka_na        ::Union{Float64,Missing}    = 13                            
    nka_k        ::Union{Float64,Missing}    =0.2                             

    
    # Part 1C. Transmembrane receptors
    syn_act     ::Union{Float64,Missing}   = missing          # [1/ms], forward rate constant                                   
    syn_deact     ::Union{Float64,Missing}   = missing         # [1/ms], backward rate constant
    syn_th ::Union{Float64,Missing} = missing

    # Part 1E. KCL cotransporter
    UKCl     ::Union{Float64,Missing} = 1f-6 # [fmol/ms/mV] the cotransporter strength 
    PWi        ::Union{Float64,Missing} = 2*1f-14   #  Water permeability

    # Part 1F. Ion related parameters
    NaCe0     ::Union{Float64,Missing} = 152   # [mMol/l = mM], Initial concentration of Sodium (o = extracellular, e=excitatory) 
    NaCi0   ::Union{Float64,Missing} = 13    # [mMol/l = mM], ICS sodium concentration of excitatory population
    KCe0      ::Union{Float64,Missing} = 3    # [mMol/l = mM], ECS potassium concentration
    KCi0    ::Union{Float64,Missing} = 145   # [mMol/l = mM], Initial ICS potassium concentration of excitatory population
    ClCe0     ::Union{Float64,Missing} = 135
    ClCi0   ::Union{Float64,Missing} = 8.0   # [mMol/l = mM], ICS chloride concentration of excitatory population
    Vi0 ::Union{Float64,Missing} = -65.5 
    

    # Part 1G. Leak permeabilities (with voltage-gated currents)
    PNaL       ::Union{Float64,Missing}  = 0.001131718785164   # [1000 mum^3/ms], Leak Na+ permeability (for -65 mV) 
    PKL        ::Union{Float64,Missing}  = 0.013554843976053   # [1000 mum^3/ms], Leak K+ permeability (for -65mV)
    PClL       ::Union{Float64,Missing}  = 0.001930766455246   # [1000 mum^3/ms], Leak Cl- permeability (for -65mV)
    

    # Part 1H. Permeabilities
    PNaG       ::Union{Float64,Missing} = 80*1f-5                               # [1000 mum^3/ms] Maximal transient Na+ permeability         
    PKG        ::Union{Float64,Missing} = 40*1f-5                               # [1000 mum^3/ms] Maximal delayed rectifier K+ permeability 
    PClG       ::Union{Float64,Missing} = 1.95*1f-5                             # [1000 mum^3/ms] Maximal voltage-gated Cl- permeability   
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
    
    X0 :: Union{Vector{Float64},Missing} = missing
    O2e_th_NKA :: Union{Float64,Missing} = 1.25
    O2e_th_vATP :: Union{Float64,Missing} = 1.5
    O2e_fac :: Union{Float64,Missing} = 0.05 # 0.1875
    min_vATP :: Union{Float64,Missing} = 0.1
    O2bath :: Union{Float64,Missing} = 2
    O2_baseline :: Union{Float64,Missing} = 1.75
    O2_alpha :: Union{Float64,Missing} = 1/6
    O2_lambda :: Union{Float64,Missing} = 1
    PvATP :: Union{Float64,Missing} = missing
    O2_diff :: Union{Float64,Missing} = 0.33*1e-3
end


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
    X0 :: Union{Array{Float64},Missing} = missing
end             


# --------------------- Helper functions ----------------------------
function MemPot(pop,x)
    NNa, NK, NCl, Vol = x
    return (pop.F/(pop.C))*(NNa+NK-NCl-pop.NAi)
end

function Nernst(pop,ext,in;z=1)
    return pop.R*pop.T/pop.F/z*log(ext/in)
end

function NKA(pop::NeuralPopSoma, NaCi, KCe, NaCe, V, O2e)
    sigmapump = 1/7*(exp(NaCe/67.3)-1) 
    fpump = 1/(1+0.1245*exp(-0.1*pop.F/pop.R/pop.T*V) +
               0.0365*sigmapump*exp(-pop.F/pop.R/pop.T*V))
    if typeof(pop).parameters[1] == Thalamus
        return fpump*(NaCi^(1.5)/(NaCi^(1.5)+pop.nka_na^(1.5)))*(KCe/(KCe+pop.nka_k))
    elseif typeof(pop).parameters[1] == Cortex 
        return 1/(1+exp((pop.O2e_th_NKA-O2e)/pop.O2e_fac))*fpump*(NaCi^(1.5)/(NaCi^(1.5)+pop.nka_na^(1.5)))*(KCe/(KCe+pop.nka_k))
    end
end

function GHK(pop,z,XCi, XCe, V)
    return (pop.F^2/pop.R/pop.T)*V*((XCi-XCe*exp(-z*pop.F*V/pop.R/pop.T))/(1-exp(-z*pop.F*V/pop.R/pop.T)));
end

function KCl(pop,KCi,KCe,ClCi,ClCe)
    (pop.R*pop.T/pop.F)*log((KCi*ClCi)/(KCe*ClCe))
end

# --------------- Somatic neural population right hand side -------------------- 
function (pop::NeuralPopSoma{A,  B})(hp::HyperParam,x,x_ECS,t,expr=nothing) where {A,B}
    NNa, NK, NCl, Wi = x
    NaCi = NNa/Wi
    KCi = NK/Wi
    ClCi = NCl/Wi
    NaCe, KCe, ClCe, We, O2e, NAe = x_ECS
    V = MemPot(pop, x)
    
    # Energy deprivation
    
    # Leaks
    INaL = pop.PNaL*GHK(pop,1,NaCi,NaCe,V)
    IKL = pop.PKL*GHK(pop,1,KCi,KCe,V)
    IClL = pop.PClL*GHK(pop,-1,ClCi,ClCe,V)

    # Gated
    alpha_n = 0.01 * (V+34.0)/( 1.0 - exp(-0.1 * (V+34.0)) ); #[no units]
    beta_n  = 0.125 * exp(-(V+44.0)/80.0);
    alpha_m = 0.1 * (V+30.0)/( 1.0 - exp(-0.1 * (V+30.0)) );
    beta_m  = 4.0 * exp(-(V+55.0)/18.0);
    alpha_h = 0.07 * exp(-(V+44.0)/20.0);
    beta_h  = 1.0/( 1.0 + exp(-0.1 * (V+14.0)) );
    m_inf = alpha_m/(alpha_m + beta_m);
    n_inf = alpha_n/(alpha_n + beta_n);
    h_inf = alpha_h/(alpha_h + beta_h);

    INaG = m_inf^3*h_inf*pop.PNaG*GHK(pop,1,NaCi,NaCe,V)
    IKG = n_inf^4*pop.PKG*GHK(pop,1,KCi,KCe,V)
    IClG = 1/(1+exp(-(V+10)/10))*pop.PClG*GHK(pop,-1,ClCi,ClCe,V)

    #Pump
    Ipump = pop.PumpStrength*NKA(pop,NaCi,KCe,NaCe,V,O2e)
    
    # KCl
    JKCl = pop.UKCl*KCl(pop,KCi,KCe,ClCi,ClCe)
    
    # Volumes
    SCi = NaCi + KCi + ClCi + pop.NAi/Wi
    SCe = NaCe + KCe + ClCe + NAe/We
    
    if expr == nothing
        return [-1/(pop.F)*(INaG + 3*Ipump + INaL);
                -1/(pop.F)*(IKG - 2*Ipump + IKL) - JKCl;
                1/(pop.F)*(IClG + IClL) - JKCl;
                pop.PWi*pop.R*pop.T*(SCi-SCe)]
    elseif expr == "NaCi"  ; return NaCi    
    elseif expr == "V"; return V
    elseif expr == "KCi"   ; return KCi                   
    elseif expr == "ClCi"  ; return ClCi                   
    elseif expr == "NaCe"  ; return NaCe                   
    elseif expr == "KCe"   ; return KCe                   
    elseif expr == "ClCe"  ; return ClCe                   
    elseif expr == "Wi"    ; return Wi                    
    elseif expr == "INaL"  ; return INaL
    elseif expr == "IKL"  ; return IKL
    elseif expr == "IClL"   ; return IClL 
    elseif expr == "alpha_n" ; return alpha_n                  
    elseif expr == "beta_n " ; return beta_n                    
    elseif expr == "alpha_m" ; return alpha_m                  
    elseif expr == "beta_m " ; return beta_m 
    elseif expr == "alpha_h" ; return alpha_h
    elseif expr == "beta_h" ; return beta_h 
    elseif expr == "m_inf" ; return m_inf
    elseif expr == "n_inf" ; return n_inf
    elseif expr == "h_inf" ; return h_inf
    elseif expr == "INaG" ; return INaG
    elseif expr == "IKG" ; return IKG
    elseif expr == "IClG" ; return IClG
    elseif expr == "Ipump" ; return Ipump
    elseif expr == "JKCl" ; return JKCl
    elseif expr == "SCi" ; return SCi
    elseif expr == "SCe" ; return SCe
    end
end

