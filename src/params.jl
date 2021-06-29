function Par(nm)
    nm.hp.ratio             = 0.8f0
    nm.Ari0             = 4*pi*(nm.rcell)^2/1000;             # [mu m^2], Membrane surface of neuron
    nm.Wi0             = (4/3)*pi*(nm.rcell)^3/1000;         # [mu m^3], Volume of spherical cell
    nm.We0             = nm.hp.ratio*(2*nm.Wi0)/(1-nm.hp.ratio)
    nm.Wtot = 2*nm.Wi0 + nm.We0
    nm.xi                = nm.F/(nm.R*nm.T);                # [1/mV]
    nm.CNa               = nm.NaCe0*nm.We0+(nm.NaCiE0+ nm.NaCiI0)*nm.Wi0;
    nm.CK                = nm.KCe0*nm.We0 +(nm.KCiE0 + nm.KCiI0)*nm.Wi0;
    nm.CCl               = nm.ClCe0*nm.We0+(nm.ClCiE0+ nm.ClCiI0)*nm.Wi0;
    nm.E_na_E_0          = 26.64f0*log(nm.NaCe0/nm.NaCiE0);       # [mV], Initial Nernst potential of sodium
    nm.E_na_I_0          = 26.64f0*log(nm.NaCe0/nm.NaCiI0);       # [mV], Initial Nernst potential of sodium of inhibitory population
    nm.E_k_E_0           = 26.64f0*log(nm.KCe0/nm.KCiE0);         # [mV], Initial Nernst potential of potassium
    nm.E_k_I_0           = 26.64f0*log(nm.KCe0/nm.KCiI0);         # [mV], Initial Nernst potential of potassium of inhibitory population

    # Get impermeants
    nm.NAi           = nm.Wi0*(-nm.C*nm.Vi0/nm.F/nm.Wi0 + nm.NaCiI0 + nm.KCiI0 - nm.ClCiI0)
    SCi              = nm.NaCiI0 + nm.KCiI0 + nm.ClCiI0 + nm.NAi/nm.Wi0
    nm.NAe           = nm.We0*(SCi -nm.NaCe0 - nm.KCe0 - nm.ClCe0)
    if nm.NAi < 0
        throw(DomainError(nm.NAi,"Impermeant NAi showuld be nonnegative"))
    elseif nm.NAe <0 
        throw(DomainError(nm.NAe,"Impermeant NAe showuld be nonnegative"))
    end

    # Get leak conductances
    # Excitatory
    X0 = [nm.NaCiE0*nm.Wi0,nm.NaCiI0*nm.Wi0,nm.KCiE0*nm.Wi0,nm.KCiI0*nm.Wi0,nm.ClCiE0*nm.Wi0,nm.ClCiI0*nm.Wi0,0.08f0,0.008f0,nm.Wi0,nm.Wi0] 
    I_TNa_E = nm(X0,0,0,"I_TNa_E")
    I_DK_E = nm(X0,0,0,"I_DK_E")
    I_GCl_E = nm(X0,0,0,"I_GCl_E")
    I_LCl_E = nm(X0,0,0,"I_LCl_E")/nm.PclL
    I_LNa_E = nm(X0,0,0,"I_LNa_E")/nm.PnaL
    I_LK_E = nm(X0,0,0,"I_LK_E")/nm.PkL
    J_KCL_E = nm(X0,0,0,"J_KCL_E")
    Ipump_E = nm(X0,0,0,"Ipump_E")
    nm.PnaL = (-I_TNa_E - 3*Ipump_E)/I_LNa_E
    nm.PkL = (-I_DK_E + 2*Ipump_E - nm.F*J_KCL_E)/I_LK_E
    nm.PclL = (-I_GCl_E + nm.F*J_KCL_E)/I_LCl_E

    if nm.PnaL < 0
        throw(DomainError(nm.PnaL,"Na leak conductance should be nonnegative"))
    elseif nm.PkL < 0
        throw(DomainError(nm.PkL,"K leak conductance should be nonnegative"))
    elseif nm.PclL < 0
        throw(DomainError(nm.PclL,"Cl leak conductance should be nonnegative"))
    end 
    
    nm.hp.excite        = [0.03f0,20.0f0,20.5f0,5f4]
    nm.X0 = X0
end

