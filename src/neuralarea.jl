Base.@kwdef mutable struct NeuralArea 
    CNa :: Union{Missing, Float64} = missing
    CK :: Union{Missing, Float64} = missing
    CCl :: Union{Missing, Float64} = missing
    We0 :: Union{Missing,Float64} = missing
    Wtot :: Union{Missing, Float64} = missing
    NAe :: Union{Missing,Float64} = missing
    pop1 :: Union{Missing, NeuralPopSoma}
    pop2 :: Union{Missing, NeuralPopSoma}
    X0 :: Union{Missing,Vector{Float64}} = missing
end 

function get_ECS(area::NeuralArea, x)
    NNa1, NK1, NCl1, Wi1, NNa2, NK2, NCl2, Wi2, O2e = x # x = [area.pop1;area.pop2]
    We = area.Wtot - Wi1 - Wi2
    NaCe = (area.CNa - NNa1 - NNa2)/We
    KCe = (area.CK - NK1 - NK2)/We
    ClCe = (area.CCl - NCl1 - NCl2)/We
    x_ECS = [NaCe; KCe; ClCe; We; O2e; area.NAe]
    return x_ECS
end

function (area::NeuralArea)(hp::HyperParam,x,t,expr=nothing)
    x_ECS = get_ECS(area, x)
    x_1 = x[1:4]
    x_2 = x[5:8]
    O2e = x[9]
    W_pop1 = x_1[4]
    W_pop2 = x_2[4]

    IvATP_1 = area.pop1.min_vATP + (1-area.pop1.min_vATP)/(1+exp((area.pop1.O2e_th-O2e)/area.pop1.O2e_fac))
    IvATP_2 = area.pop2.min_vATP + (1-area.pop2.min_vATP)/(1+exp((area.pop2.O2e_th-O2e)/area.pop2.O2e_fac))

    Ipump_1 = area.pop1.PumpStrength*NKA(area.pop1,x_1[1]/x_1[4],x_ECS[2],x_ECS[1],MemPot(area.pop1,x_1),O2e)
    Ipump_2 = area.pop2.PumpStrength*NKA(area.pop2,x_2[1]/x_2[4],x_ECS[2],x_ECS[1],MemPot(area.pop2,x_2),O2e)

    
    IBlock = 1/(1+exp(hp.beta1*(t-hp.tstart))) + 1/(1+exp(-hp.beta2*(t-hp.tend)))
    IBlock = hp.perc + (1-hp.perc)*IBlock
    O2bath = IBlock*area.pop1.O2bath

    O2_rhs = -area.pop1.O2_alpha*area.pop1.O2_lambda*(1/area.pop1.F)*(Ipump_1/W_pop1 + Ipump_2/W_pop2 + area.pop1.PvATP*IvATP_1/W_pop1 + area.pop2.PvATP*IvATP_2/W_pop2) + area.pop1.O2_diff*(O2bath - O2e)

    if expr == nothing
        return [area.pop1(hp,x_1,x_ECS,t); area.pop2(hp,x_2,x_ECS,t); O2_rhs]
    elseif expr[end] == 'E'
        if typeof(area.pop1).parameters[2] == Excitatory
            return area.pop1(hp,x_1,x_ECS,t,expr[1:end-1])
        else
            return area.pop2(hp,x_2,x_ECS,t,expr[1:end-1])
        end
    elseif expr[end] == 'I'
        if typeof(area.pop1).parameters[2] == Inhibitory
            return area.pop1(hp,x_1,x_ECS,t,expr[1:end-1])
        else
            return area.pop2(hp,x_2,x_ECS,t,expr[1:end-1])
        end
    end 
end

