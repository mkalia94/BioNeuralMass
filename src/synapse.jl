Base.@kwdef mutable struct BioNM
    conn :: Union{Matrix,Missing} = missing  # CAREFUL: size(conn) needs to be 2*(n_areas) x 2*(n_areas)
    sol :: Union{Missing,ODESolution} = missing
    hp :: Union{HyperParam,Missing} = missing
    areas :: Union{Vector{NeuralArea},Missing} = missing
end

function syn_current_pop(pop::NeuralPopSoma,x,x_ECS,t,type="Excitatory")
    if type == "Excitatory"
        NaCe = x_ECS[1]
        NaCi = x[1]/x[4]
        return -(MemPot(pop,x)-Nernst(pop,NaCe,NaCi))
        #return (NaCe-NaCi)
    else
        ClCe = x_ECS[3]
        ClCi = x[3]/x[4]
        #return -(ClCe-ClCi)
        return -(MemPot(pop,x)-Nernst(pop,ClCe,ClCi,z=-1))
    end
end

function syn_current(area::NeuralArea,x,t)
    x_ECS = get_ECS(area, x)
    syn11 = syn_current_pop(area.pop1,x[1:4],x_ECS,t,"Excitatory")
    syn12 = syn_current_pop(area.pop1,x[1:4],x_ECS,t,"Inhibitory")
    syn21 = syn_current_pop(area.pop2,x[5:8],x_ECS,t,"Excitatory")
    syn22 = syn_current_pop(area.pop2,x[5:8],x_ECS,t,"Inhibitory")
    return [syn11 syn12; syn21 syn22]
end

function syn_act(pop::NeuralPopSoma,FR)
    return pop.syn_act*(FR)/(FR + pop.syn_th);
end

# Firing Rate
function firing_rate(pop::NeuralPopSoma,I,EK,ENa)
    Ith     = pop.p00 + pop.p10*EK + pop.p01*ENa + pop.p20*(EK)^2 + pop.p11*EK*ENa + pop.p02*(ENa)^2 + pop.p30*(EK)^3 + pop.p21*(EK)^2*ENa + pop.p12*EK*(ENa)^2 + pop.p03*(ENa)^3;  # 1st threshold
    kappa   = pop.q00 + pop.q10*EK + pop.q01*ENa + pop.q20*(EK)^2 + pop.q11*EK*ENa + pop.q02*(ENa)^2 + pop.q30*(EK)^3 + pop.q21*(EK)^2*ENa + pop.q12*EK*(ENa)^2 + pop.q03*(ENa)^3;  #
    Ith2    = pop.r00 + pop.r10*EK + pop.r01*ENa + pop.r20*(EK)^2 + pop.r11*EK*ENa + pop.r02*(ENa)^2 + pop.r30*(EK)^3 + pop.r21*(EK)^2*ENa + pop.r12*EK*(ENa)^2 + pop.r03*(ENa)^3;  # 2nd threshold
    fun = (x)->((1/(pop.sigma*sqrt(2*pi)))*exp(-(I-x)^2/(2*pop.sigma^2)))*(kappa*sqrt(max(0,x-Ith)))*(1-(1+sign(x-Ith2))/2);
    #prob1 = QuadratureProblem(fun,0,Inf)
    FR,_ = quadgk(fun,-Inf,Inf); # NEEDS: QuadGK
    #FR = solve(prob1, QuadGKJL(),reltol=1e-3,abstol=1e-3).u 
    # FR = (kappa*sqrt(max(0,I-Ith)))*(1-(1+sign(I-Ith2))/2);
   return FR
end

function get_EEG(areas::Vector{NeuralArea},syn_curr)
    ctr = 1
    for area in areas
        if typeof(area.pop1).parameters[1] == Cortex && typeof(area.pop1).parameters[2] == Excitatory
            return syn_curr[(ctr-1)*2+1]
        elseif typeof(area.pop2).parameters[1] == Cortex && typeof(area.pop2).parameters[2] == Excitatory
            return syn_curr[ctr*2]
        end
        ctr = ctr+1
    end
end


function syn_var_RHS_pop(hp::HyperParam,pop::NeuralPopSoma,rsyn,x,x_ECS,t,I_Ext)
    NaCe = x_ECS[1]
    KCe = x_ECS[2]
    O2e = x_ECS[5]
    NaCi = x[1]/x[4]
    KCi = x[2]/x[4]
    ENa = Nernst(pop,NaCe,NaCi)
    EK = Nernst(pop,KCe,KCi)
    I = I_Ext - pop(hp,x,x_ECS,t,"Ipump")
    FR = firing_rate(pop, I,EK, ENa)
    rsyn_act = syn_act(pop,FR)
    block = hp.min_vATP + (1-hp.min_vATP)/(1+exp((hp.O2e_th_vATP-O2e)/pop.O2e_fac))
    if typeof(pop).parameters[1] == Thalamus
        return rsyn_act*(1-rsyn) - pop.syn_deact*rsyn
    else
        return block*(rsyn_act*(1-rsyn) - pop.syn_deact*rsyn)
    end
    end
    

function syn_var_RHS(hp::HyperParam,area::NeuralArea,rsyn,x,t,I_Ext)
    x_ECS = get_ECS(area,x)
    syn1 = syn_var_RHS_pop(hp,area.pop1,rsyn[1],x[1:4],x_ECS,t,I_Ext[1])
    syn2 = syn_var_RHS_pop(hp,area.pop2,rsyn[2],x[5:8],x_ECS,t,I_Ext[2])
    return [syn1;syn2]
end

function (nm::BioNM)(dx,x,p,t,expr=nothing)
    
    # Note: x_Bio = [NNa,NK,NCl,W], x_syn = [rsyn]; x_area = [x_Bio_1,x_Bio_2,x_syn_1,x_syn_2]; x = [x_area_1...] 

    # External stimulation
    if !ismissing(nm.hp.excite)
        dur_ = nm.hp.excite[3]-nm.hp.excite[2]
        curr_ = nm.hp.excite[1]
        fac_ = nm.hp.excite[4]
        IExcite = curr_*(squarewave1((t-nm.hp.excite[2])/fac_,dur_/fac_) + 1)/2 
    else
        IExcite = 0
    end

    # Handle synapses
    rsyn = vcat([x[11*i-1:11*i] for i in 1:div(length(x),11)]...)
    # syn_curr = zeros(2*length(nm.areas),2)
    
    syn_curr = repeat(vcat([syn_current(nm.areas[ctr],x[(ctr-1)*11+1:ctr*11],t) for ctr in 1:length(nm.areas)]...),1,length(nm.areas))
    
    # ctr = 1
    # for area in nm.areas
    #     x_1 = x[(ctr-1)*10+1:ctr*10]
    #     syn_curr[(ctr-1)*2 + 1:ctr*2,:] = syn_current(area,x_1,t)
    #     ctr = ctr + 1
    # end
    
    # for i in 1:(length(nm.areas)-1)
    #     syn_curr = hcat(syn_curr,syn_curr)
    # end
    
    syn_curr_full = rsyn' .* (nm .* syn_curr)
    syn_curr_full = [sum(syn_curr_full[1:2:end]) sum(syn_curr_full[2:2:end])] # [Excitatory Inhibitory]
    syn_curr = (nm.conn .* syn_curr)*rsyn
    
    # Finally, generate RHS
    ctr = 1
    I_EEG = nothing
    FR = zeros(4)

    for area in nm.areas
        x_1 = x[(ctr-1)*11+1:ctr*11]
        rsyn_1 = rsyn[(ctr-1)*2+1:ctr*2] 
        syn_curr_1 = syn_curr[(ctr-1)*2+1:ctr*2]
        syn_curr_full_1 = syn_curr_full[(ctr-1)*2+1:ctr*2,:]
        
        if ctr == 1 && typeof(area.pop1).parameters[2] == Excitatory
            I_Ext = [IExcite;0] .+ syn_curr_1
        elseif ctr == 1 && typeof(area.pop2).parameters[2] == Excitatory
            I_Ext = [0;IExcite] .+ syn_curr_1
        else
            I_Ext = [0;0] .+ syn_curr_1
        end
        
        if !nm.hp.synapseoff && expr == nothing
            dx[(ctr-1)*11+1:ctr*11] = [area(nm.hp,x_1,syn_curr_full_1,t); syn_var_RHS(nm.hp,area,rsyn_1,x_1,t,I_Ext)]
        elseif expr == "FiringRate"
            x_ECS = get_ECS(area,x_1)
            EK1 = area.pop1.R*area.pop1.T/area.pop1.F*(log(x_ECS[2]/x_1[2]*x_1[4]))
            ENa1 = area.pop1.R*area.pop1.T/area.pop1.F*(log(x_ECS[1]/x_1[1]*x_1[4]))
            EK2 = area.pop1.R*area.pop1.T/area.pop1.F*(log(x_ECS[2]/x_1[6]*x_1[8]))
            ENa2 = area.pop1.R*area.pop1.T/area.pop1.F*(log(x_ECS[2]/x_1[5]*x_1[8]))
            Ipump = [area.pop1(nm.hp,x_1[1:4],x_ECS,t,"Ipump");
                     area.pop2(nm.hp,x_1[5:8],x_ECS,t,"Ipump")]
            II = I_Ext .- Ipump 
            FR[(ctr-1)*2+1:ctr*2] = [firing_rate(area.pop1,II[1],EK1,ENa1);
                                     firing_rate(area.pop2,II[2],EK2,ENa2)]
        else
            dx[(ctr-1)*11+1:ctr*11] = [area(nm.hp,x_1,t); 0; 0]
        end

        ctr = ctr+1
    end
    

    if expr == "EEGraw"; return get_EEG(nm.areas,syn_curr)
    elseif expr == "rsyn"; return rsyn
    elseif expr == "FiringRate"; return FR
    elseif expr == "syn_curr"; return syn_curr
    elseif expr == "IExcite"; return IExcite
    else; return dx
    end
   

end



function (nm::BioNM)(expr::String=nothing)
    if ismissing(nm.sol)
        throw(TypeError(nm.sol),ODESolution,missing)
    else
        tmp_ = hcat([nm(zeros(size(nm.sol.u[i])),nm.sol.u[i],0.0,nm.sol.t[i],expr) for i in 1:length(nm.sol.t)]...)
        return tmp_'
    end
end

function (nm::BioNM)(area::Symbol,expr::String=nothing)
    ctr = 1
    for ar_ in nm.areas
        if typeof(ar_.pop1).parameters[1] == eval(area)
            if ismissing(nm.sol)
                throw(TypeError(nm.sol),ODESolution,missing)
            else
                tmp_ =  hcat([ar_(nm.hp,nm.sol.u[i][(ctr-1)*11+1:11*ctr],nm.sol.t[i],expr) for i in 1:length(nm.sol.t)]...)
                # return reshape(tmp_,size(tmp_)[2],size(tmp_)[1])
                return tmp_'
            end
        end
        ctr = ctr + 1
    end
end




