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
    NNa1, NK1, NCl1, Wi1, NNa2, NK2, NCl2, Wi2 = x # x = [area.pop1;area.pop2]
    We = area.Wtot - Wi1 - Wi2
    NaCe = (area.CNa - NNa1 - NNa2)/We
    KCe = (area.CK - NK1 - NK2)/We
    ClCe = (area.CCl - NCl1 - NCl2)/We
    x_ECS = [NaCe; KCe; ClCe; We; area.NAe]
    return x_ECS
end

function (area::NeuralArea)(hp::HyperParam,x,t,expr=nothing)
    x_ECS = get_ECS(area, x)
    x_1 = x[1:4]
    x_2 = x[5:8]
    if expr == nothing
        return [area.pop1(hp,x_1,x_ECS,t); area.pop2(hp,x_2,x_ECS,t)]
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

