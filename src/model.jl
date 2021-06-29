# AMPA
function ampa(nm,FR_E,FR_ext)
    FR = FR_E + FR_ext;
    a_ampa = nm.a_ampa*(FR)/(FR + 15);
end

#GABA
function gaba(nm,FR_I,FR_ext)
    FR = FR_I + FR_ext;
    a_gaba = nm.a_gaba*(FR)/(FR + 100);
end

# Firing Rate
function Firing_Rate(nm,I,Ek,Ena)
    Ith     = nm.p00 + nm.p10*Ek + nm.p01*Ena + nm.p20*(Ek)^2 + nm.p11*Ek*Ena + nm.p02*(Ena)^2 + nm.p30*(Ek)^3 + nm.p21*(Ek)^2*Ena + nm.p12*Ek*(Ena)^2 + nm.p03*(Ena)^3;  # 1st threshold
    kappa   = nm.q00 + nm.q10*Ek + nm.q01*Ena + nm.q20*(Ek)^2 + nm.q11*Ek*Ena + nm.q02*(Ena)^2 + nm.q30*(Ek)^3 + nm.q21*(Ek)^2*Ena + nm.q12*Ek*(Ena)^2 + nm.q03*(Ena)^3;  #
    Ith2    = nm.r00 + nm.r10*Ek + nm.r01*Ena + nm.r20*(Ek)^2 + nm.r11*Ek*Ena + nm.r02*(Ena)^2 + nm.r30*(Ek)^3 + nm.r21*(Ek)^2*Ena + nm.r12*Ek*(Ena)^2 + nm.r03*(Ena)^3;  # 2nd threshold
    fun = x->((1/(nm.sigma*sqrt(2f0*pi)))*exp(-(I-x)^2/(nm.sigma)))*(kappa*sqrt(x-Ith))*(1-(1+sign(x-Ith2))/2);
    FR,_ = quadgk(fun,Ith,300f0); # NEEDS: QuadGK
    return FR
end

function (nm::BioNM)(x,p,t,expr=nothing)
    
    NNaE = x[1]
    NNaI = x[2]
    NKE = x[3]
    NKI = x[4]
    NClE = x[5]
    NClI = x[6]
    rAMPA = x[7]
    rGABA = x[8]
    WiE = x[9]
    WiI = x[10]
    
    NaCiI = NNaI/WiI
    NaCiE = NNaE/WiE
    KCiI = NKI/WiI
    KCiE = NKE/WiE
    ClCiI = NClI/WiI
    ClCiE = NClE/WiE

    We = nm.Wtot - WiI - WiE

    NaCe = (nm.CNa - NNaI - NNaE)/We
    KCe = (nm.CK - NKI - NKE)/We
    ClCe = (nm.CCl - NClI - NClE)/We

    IBlock = 1/(1+exp(nm.hp.beta1*(t-nm.hp.tstart))) + 1/(1+exp(-nm.hp.beta2*(t-nm.hp.tend)))
    IBlock = nm.hp.perc + (1-nm.hp.perc)*IBlock
    
    if !ismissing(nm.hp.excite)
        dur_ = nm.hp.excite[3]-nm.hp.excite[2]
        curr_ = nm.hp.excite[1]
        fac_ = nm.hp.excite[4]
        IExcite = (squarewave1((t-nm.hp.excite[2])/fac_,dur_/fac_) + 1f0)/2f0 |> Float32
    else
        IExcite = 0f0
    end 


    # Voltage
    V_E   = (nm.F/(nm.C))*(WiE)*(NaCiE+KCiE-ClCiE-nm.NAi/WiE);
    V_I   = (nm.F/(nm.C))*(WiI)*(NaCiI+KCiI-ClCiI-nm.NAi/WiI);

    # Nernst potentials (Excitatory population)
    E_na_E                = 26.64f0*log(NaCe/NaCiE);          # [mV], Nernst potential sodium
    E_k_E                 = 26.64f0*log(KCe/KCiE);           # [mV], Nernst potential potassium
    E_cl_E                = 26.64f0*log(ClCiE/ClCe);          # [mV], Nernst potential chloride

    # Nernst potentials (Inhibitory population)
    E_na_I                = 26.64f0*log(NaCe/NaCiI);
    E_k_I                 = 26.64f0*log(KCe/KCiI);
    E_cl_I                = 26.64f0*log(ClCiI/ClCe);

    # Equilibrium potentials (variable)
    #nm.V_ampa_eq_E           = (E_na_E + E_k_E)/2;           # [mV], equilibrium potential (AMPA)
    #nm.V_gaba_eq_E           =  E_cl_E;                      # [mV], equilibrium potential (GABA A)
    #nm.V_ampa_eq_I           = (E_na_I + E_k_I)/2;           # [mV], equilibrium potential (AMPA)
    #nm.V_gaba_eq_I           =  E_cl_I;                      # [mV], equilibrium potential (GABA A)
    I_E = (Float32(1/500)*(nm.F*nm.Wi0)/(nm.C*nm.Ari0))*((nm.g_ampa_E)*rAMPA*(NaCe - NaCiE) - nm.g_gaba_E*rGABA*(ClCe - ClCiE))
    I_I =  (Float32(1/500)*(nm.F*nm.Wi0)/(nm.C*nm.Ari0))*(nm.g_ampa_I*rAMPA*(NaCe - NaCiI) - nm.g_gaba_I*rGABA*(ClCe - ClCiI))
    I_E = (nm.g_ampa_E*rAMPA*(NaCe - NaCiE) - nm.g_gaba_E*rGABA*(ClCe - ClCiE))
    I_I = (nm.g_ampa_I*rAMPA*(NaCe - NaCiI) - nm.g_gaba_I*rGABA*(ClCe - ClCiI))


   # Ion pump
    sigmapump = 1/7*(exp(NaCe/67.3)-1) |> Float32
    fpump_I = 1/(1+0.1245*exp(-0.1*nm.F/nm.R/nm.T*V_I) +
               0.0365*sigmapump*exp(-nm.F/nm.R/nm.T*V_I)) |> Float32
    fpump_E = 1/(1+0.1245*exp(-0.1*nm.F/nm.R/nm.T*V_E) +
               0.0365*sigmapump*exp(-nm.F/nm.R/nm.T*V_E)) |> Float32
    Ipump_E = IBlock*fpump_E*nm.rho_p_E*(
                                       NaCiE^(1.5)/(NaCiE^(1.5)+nm.nka_na^(1.5)))*(KCe/(KCe+nm.nka_k))
    Ipump_I = IBlock*fpump_I*nm.rho_p_E*(
                                       NaCiI^(1.5)/(NaCiI^(1.5)+nm.nka_na^(1.5)))*(KCe/(KCe+nm.nka_k))
    FR_E                  = Firing_Rate(nm,I_E-Ipump_E,E_k_E,E_na_E);                                     # [1/ms], Firing Rate
    FR_I                  = Firing_Rate(nm,I_I-Ipump_I,E_k_I,E_na_I);                                     # [1/ms], Firing Rate

    # Leak currents (excitatory population)
    I_LNa_E               = nm.PnaL*(nm.F*nm.xi)*V_E*((NaCiE-NaCe*exp(-nm.xi*V_E))/(1-exp(-nm.xi*V_E)));  # Leak current sodium
    I_LK_E                = nm.PkL*(nm.F*nm.xi)*V_E*((KCiE-KCe*exp(-nm.xi*V_E))/(1-exp(-nm.xi*V_E)));    # Leak current potassium
    I_LCl_E               = nm.PclL*(nm.F*nm.xi)*V_E*((ClCiE-ClCe*exp(nm.xi*V_E))/(1-exp(nm.xi*V_E)));    # Leak current Chloride z = -1

    # Leak currents (inhinitatory population)
    I_LNa_I               = nm.PnaL*(nm.F*nm.xi)*V_I*((NaCiI-NaCe*exp(-nm.xi*V_I))/(1-exp(-nm.xi*V_I)));  # Leak current sodium
    I_LK_I                = nm.PkL*(nm.F*nm.xi)*V_I*((KCiI-KCe*exp(-nm.xi*V_I))/(1-exp(-nm.xi*V_I)));    # Leak current potassium
    I_LCl_I               = nm.PclL*(nm.F*nm.xi)*V_I*((ClCiI-ClCe*exp(nm.xi*V_I))/(1-exp(nm.xi*V_I)));     # Leak current Chloride z = -1

    # Voltage-gated sodium/potassium/chloride current (excitatory population)
    alpha_n     = 0.01f0 * (V_E+34.0f0)/( 1.0f0 - exp(-0.1f0 * (V_E+34.0f0)) ); #[no units]
    beta_n      = 0.125f0 * exp(-(V_E+44.0f0)/80.0f0);
    alpha_m     = 0.1f0 * (V_E+30.0f0)/( 1.0f0 - exp(-0.1f0 * (V_E+30.0f0)) );
    beta_m      = 4.0f0 * exp(-(V_E+55.0f0)/18.0f0);
    alpha_h     = 0.07f0 * exp(-(V_E+44.0f0)/20.0f0);
    beta_h      = 1.0f0/( 1.0f0 + exp(-0.1f0 * (V_E+14.0f0)) );

    m_inf       = alpha_m/(alpha_m + beta_m);
    n_inf       = alpha_n/(alpha_n + beta_n);
    h_inf       = alpha_h/(alpha_h + beta_h);

    I_TNa_E     = nm.PnaT*m_inf^3*h_inf*(nm.F*nm.xi)*V_E*((NaCiE-NaCe*exp(-nm.xi*V_E))/(1-exp(-nm.xi*V_E)));            # Voltage gated sodium current
    I_DK_E      = nm.PkD*n_inf^4*(nm.F*nm.xi)*V_E*((KCiE-KCe*exp(-nm.xi*V_E))/(1-exp(-nm.xi*V_E)));                    # Voltage gated potassium current
    I_GCl_E     = (nm.PclG/(1+exp(-(V_E+10f0)/10f0)))*(nm.F*nm.xi)*V_E*((ClCiE-ClCe*exp(nm.xi*V_E))/(1-exp(nm.xi*V_E)));    # Voltage gated chloride current

    # Voltage-gated sodium/potassium/chloride current (inhibitatory population)
    alpha_n     = 0.01f0 * (V_I+34.0f0)/( 1.0f0 - exp(-0.1f0 * (V_I+34.0f0)) ); #[no units]
    beta_n      = 0.125f0 * exp(-(V_I+44.0f0)/80.0f0);
    alpha_m     = 0.1f0 * (V_I+30.0f0)/( 1.0f0 - exp(-0.1f0 * (V_I+30.0f0)) );
    beta_m      = 4.0f0 * exp(-(V_I+55.0f0)/18.0f0);
    alpha_h     = 0.07f0 * exp(-(V_I+44.0f0)/20.0f0);
    beta_h      = 1.0f0/( 1.0f0 + exp(-0.1f0 * (V_I+14.0f0)) );
    m_inf       = alpha_m/(alpha_m + beta_m);
    n_inf       = alpha_n/(alpha_n + beta_n);
    h_inf       = alpha_h/(alpha_h + beta_h);

    I_TNa_I               = nm.PnaT*m_inf^3*h_inf*(nm.F*nm.xi)*V_I*((NaCiI-NaCe*exp(-nm.xi*V_I))/(1-exp(-nm.xi*V_I)));    # Voltage gated sodium current
    I_DK_I                = nm.PkD*n_inf^4*(nm.F*nm.xi)*V_I*((KCiI-KCe*exp(-nm.xi*V_I))/(1-exp(-nm.xi*V_I)));  # Voltage gated potassium current
    I_GCl_I               = (nm.PclG/(1+exp(-(V_I+10f0)/10f0)))*(nm.F*nm.xi)*V_I*((ClCiI-ClCe*exp(nm.xi*V_I))/(1-exp(nm.xi*V_I)));

    # KLC costransporter
    J_KCL_E   = nm.U_kcl*(1/nm.xi)*log((KCiE*ClCiE)/(KCe*ClCe));
    J_KCL_I   = nm.U_kcl*(1/nm.xi)*log((KCiI*ClCiI)/(KCe*ClCe));

    # Recovery (Block voltage gated sodium channels)
    # if t > nm.Trecovery
    #     I_TNa_E = 0;
    #     I_TNa_I = 0;
    # end
    
    # Water permeability
    SCiI = NaCiI + KCiI + ClCiI + nm.NAi/WiI
    SCiE = NaCiE + KCiE + ClCiE + nm.NAi/WiE
    SCe = NaCe + KCe + ClCe + nm.NAe/We
   

    # State vector

dX = [(1/(nm.F))*(-I_LNa_E-I_TNa_E - 3*Ipump_E)  + (1/(nm.F))*(nm.g_ampa_E*rAMPA*(NaCe - NaCiE));          # Na_i excitatory
           (1/(nm.F))*(-I_LNa_I-I_TNa_I - 3*Ipump_I)  + (1/(nm.F))*(nm.g_ampa_I*rAMPA*(NaCe - NaCiI));          # Na_i inhibitory
           (1/(nm.F))*(-I_LK_E-I_DK_E   + 2*Ipump_E)    - J_KCL_E;                         # dK_i/dt excitatory
           (1/(nm.F))*(-I_LK_I-I_DK_I   + 2*Ipump_I)    - J_KCL_I;                         # dK_i/dt inhibitory
           (1/(nm.F))*(I_LCl_E+I_GCl_E)   - J_KCL_E + (1/(nm.F*nm.C))*(nm.g_gaba_E*rGABA*(ClCe - ClCiE));                        # dCl_i/dt excitatory
           (1/(nm.F))*(I_LCl_I+I_GCl_I)   - J_KCL_I + (1/(nm.F*nm.C))*(nm.g_gaba_I*rGABA*(ClCe - ClCiI));                   # dK_i/dt  inhibitory
           ampa(nm,FR_E,IExcite)*(1-rAMPA)-nm.b_ampa*rAMPA;                                                                                 # dr_AMPA/dt
           gaba(nm,FR_I,0)*(1-rGABA)-nm.b_gaba*rGABA;
           nm.PWi*nm.R*nm.T*(SCiE-SCe);
           nm.PWi*nm.R*nm.T*(SCiI-SCe)]                                                                                  # dr_GABA/dt
    
 dX = dX*60f0*1f3
    
    if expr == nothing      ; return  dX
    elseif expr == "NaCiE"  ; return  NaCiE                   
    elseif expr == "NaCiI"  ; return  NaCiI    
    elseif expr == "KCiE"   ; return  KCiE     
    elseif expr == "KCiI"   ; return  KCiI     
    elseif expr == "ClCiE"  ; return  ClCiE    
    elseif expr == "ClCiI"  ; return  ClCiI    
    elseif expr == "rAMPA"  ; return  rAMPA    
    elseif expr == "rGABA"  ; return  rGABA    
    elseif expr == "IBlock" ; return  IBlock   
    elseif expr == "IBlock" ; return  IBlock     
    elseif expr == "dur_"   ; return  dur_     
    elseif expr == "curr_"  ; return  curr_    
    elseif expr == "fac_"   ; return  fac_     
    elseif expr == "IExcite"; return  IExcite      
    elseif expr == "NaCe"   ; return  NaCe          
    elseif expr == "KCe"    ; return  KCe           
    elseif expr == "ClCe"   ; return  ClCe          
    elseif expr == "V_E"    ; return  V_E      
    elseif expr == "V_I"    ; return  V_I      
    elseif expr == "E_na_E" ; return  E_na_E        
    elseif expr == "E_k_E"  ; return  E_k_E         
    elseif expr == "E_cl_E" ; return  E_cl_E        
    elseif expr == "E_na_I" ; return  E_na_I        
    elseif expr == "E_k_I"  ; return  E_k_I         
    elseif expr == "E_cl_I" ; return  E_cl_I   
    elseif expr == "FR_E"   ; return  FR_E        
    elseif expr == "FR_I"   ; return  FR_I        
    elseif expr == "Ipump_E"; return  Ipump_E     
    elseif expr == "Ipump_I"; return  Ipump_I     
    elseif expr == "I_LNa_E"; return  I_LNa_E     
    elseif expr == "I_LK_E" ; return  I_LK_E      
    elseif expr == "I_LCl_E"; return  I_LCl_E     
    elseif expr == "I_LNa_I"; return  I_LNa_I     
    elseif expr == "I_LK_I" ; return  I_LK_I      
    elseif expr == "I_LCl_I"; return  I_LCl_I     
    elseif expr == "alpha_n"; return  alpha_n     
    elseif expr == "beta_n" ; return  beta_n      
    elseif expr == "alpha_m"; return  alpha_m     
    elseif expr == "beta_m" ; return  beta_m      
    elseif expr == "alpha_h"; return  alpha_h     
    elseif expr == "beta_h" ; return  beta_h      
    elseif expr == "m_inf"  ; return  m_inf       
    elseif expr == "n_inf"  ; return  n_inf       
    elseif expr == "h_inf"  ; return  h_inf       
    elseif expr == "I_TNa_E"; return  I_TNa_E     
    elseif expr == "I_DK_E" ; return  I_DK_E      
    elseif expr == "I_GCl_E"; return  I_GCl_E     
    elseif expr == "alpha_n"; return  alpha_n     
    elseif expr == "beta_n" ; return  beta_n      
    elseif expr == "alpha_m"; return  alpha_m     
    elseif expr == "beta_m" ; return  beta_m      
    elseif expr == "alpha_h"; return  alpha_h     
    elseif expr == "beta_h" ; return  beta_h      
    elseif expr == "m_inf"  ; return  m_inf       
    elseif expr == "n_inf"  ; return  n_inf       
    elseif expr == "h_inf"  ; return  h_inf       
    elseif expr == "I_TNa_I"; return  I_TNa_I     
    elseif expr == "I_DK_I" ; return  I_DK_I      
    elseif expr == "I_GCl_I"; return  I_GCl_I     
    elseif expr == "J_KCL_E"; return  J_KCL_E     
    elseif expr == "J_KCL_I"; return  J_KCL_I
    elseif expr == "WiI"    ; return  WiI
    elseif expr == "WiE"    ; return  WiE
    elseif expr == "We"    ; return  We
    elseif expr == "I_E"    ; return I_E
    elseif expr == "I_I"     ; return I_I
end
end

function (nm::BioNM)(expr::String=nothing)
    if ismissing(nm.sol)
        throw(TypeError(nm.sol),ODESolution,missing)
    else
        return nm.(nm.sol.u,0.0,nm.sol.t,expr)
    end
end

