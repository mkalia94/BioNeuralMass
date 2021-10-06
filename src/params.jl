function NeuralArea(hp::HyperParam,pop1::NeuralPopSoma,pop2::NeuralPopSoma)
     
    # Get conservation constants
    We0             = hp.ratio*(pop1.Wi0 + pop2.Wi0)/(1-hp.ratio)
    Wtot            = pop1.Wi0 + pop2.Wi0 + We0
    CNa             = pop1.NaCe0*We0+pop1.NaCi0*pop1.Wi0 + pop2.NaCi0*pop2.Wi0;
    CK              = pop1.KCe0*We0 +pop1.KCi0*pop1.Wi0+ pop2.KCi0*pop2.Wi0;
    CCl             = pop1.ClCe0*We0+pop1.ClCi0*pop1.Wi0+ pop2.ClCi0*pop2.Wi0;
    
    println("Neural Area successfully generated...")
    NeuralArea(Wtot = Wtot,We0 = We0, CNa = CNa, CK = CK, CCl = CCl, pop1 = pop1, pop2 = pop2)    
end 

function Par(hp::HyperParam,area::NeuralArea,pop::NeuralPopSoma)

    # Get impermeants
    pop.NAi           = pop.Wi0*(-pop.C*pop.Vi0/pop.F/pop.Wi0 + pop.NaCi0 + pop.KCi0 - pop.ClCi0)
    SCi              = pop.NaCi0 + pop.KCi0 + pop.ClCi0 + pop.NAi/pop.Wi0
    area.NAe           = area.We0*(SCi -pop.NaCe0 - pop.KCe0 - pop.ClCe0)
    if pop.NAi < 0
        throw(DomainError(pop.NAi,"Impermeant NAi showuld be nonnegative"))
    elseif area.NAe <0 
        throw(DomainError(area.NAe,"Impermeant area.NAe showuld be nonnegative"))
    end
    
    # Pop 1
    X0 = [pop.NaCi0*pop.Wi0,pop.KCi0*pop.Wi0,pop.ClCi0*pop.Wi0,pop.Wi0] 
    X_ECS = [pop.NaCe0, pop.KCe0,pop.ClCe0, area.We0, pop.O2_baseline, area.NAe]
    INaG = pop(hp,X0,X_ECS, 0,"INaG")
    IKG = pop(hp,X0,X_ECS,0,"IKG")
    IClG = pop(hp,X0,X_ECS,0,"IClG")
    IClL = pop(hp,X0,X_ECS,0,"IClL")/pop.PClL
    INaL = pop(hp,X0,X_ECS,0,"INaL")/pop.PNaL
    IKL = pop(hp,X0,X_ECS,0,"IKL")/pop.PKL
    JKCl = pop(hp,X0,X_ECS,0,"JKCl")
    Ipump = pop(hp,X0,X_ECS,0,"Ipump")
    pop.PNaL = (-INaG - 3*Ipump)/INaL
    pop.PKL = (-IKG + 2*Ipump - pop.F*JKCl)/IKL
    pop.PClL = (-IClG + pop.F*JKCl)/IClL

    IvATP = pop.min_vATP + (1-pop.min_vATP)/(1+exp((pop.O2e_th_vATP-pop.O2_baseline)/pop.O2e_fac))
    pop.PvATP = Ipump/(70*IvATP)
    pop.O2_diff = (pop.O2_alpha*pop.O2_lambda*(1/pop.F)*(2*Ipump/pop.Wi0+ 2*pop.PvATP*IvATP/pop.Wi0))/(pop.O2bath-pop.O2_baseline) 


    if pop.PNaL < 0
        throw(DomainError(pop.PNaL,"Na leak conductance should be nonnegative in Pop 1"))
    elseif pop.PKL < 0
        throw(DomainError(pop.PKL,"K leak conductance should be nonnegative in Pop 1"))
    elseif pop.PClL < 0
        throw(DomainError(pop.PClL,"Cl leak conductance should be nonnegative in Pop 1"))
    end

    if typeof(pop).parameters[2] == Excitatory
        pop.syn_act = 12.5*2
        pop.syn_deact = 3.0*2
        pop.syn_th = 0.2
    else
        pop.syn_act = 5.0*2  
        pop.syn_deact = 0.03*2
        pop.syn_th = 0.5
    end

    pop.X0 = X0

end

function Par(hp::HyperParam,area::NeuralArea)
    Par(hp,area,area.pop1)
    Par(hp,area,area.pop2)
    area.X0 = [area.pop1.X0; area.pop2.X0; area.pop1.O2_baseline; 0.008; 0.008]
end
