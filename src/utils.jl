import DifferentialEquations: solve
import Plots: plot
using SignalAnalysis

function solve(nm::BioNM,alg=Tsit5();kw...)
    X0 = []
    for area in nm.areas
        push!(X0,area.X0)
    end
    nm.hp.X0 = vcat(X0...)
    nm.hp.X0 = nm.hp.X0 + 1e-7*rand(length(nm.hp.X0))
    prob = ODEProblem((dx,x,p,t)->nm(dx,x,p,t),nm.hp.X0,[0.0f0,nm.hp.tfinal]) 
    sol = solve(prob,alg;kw...)
    nm.sol = sol
end

function plot(nm::BioNM,name=nothing;kw...)
    ctr = 1
    for area in nm.areas
        t = nm.sol.t
        p = []
        area_type = Symbol(typeof(area.pop1).parameters[1])

        p1 = plot(t,nm(area_type,"NaCiI"),ylabel=latexstring("[\\textrm{Na}^+]"),label="",lc=:blue)
        plot!(t,nm(area_type,"NaCiE"),label="",lc=:red)
        push!(p,p1)

        p2 = plot(t,nm(area_type,"KCiI"),ylabel=latexstring("[\\textrm{K}^+]"),label="",lc=:blue)
        plot!(t,nm(area_type,"KCiE"),label="",lc=:red)
        push!(p,p2)

        p3 = plot(t,nm(area_type,"ClCiI"),ylabel=latexstring("[\\textrm{Cl}^-]"),label="",lc=:blue)
        plot!(t,nm(area_type,"ClCiE"),label="",lc=:red)
        push!(p,p3)

        p4 = plot(t,nm(area_type,"WiI"),ylabel=latexstring("\\textrm{Vol.}"),label="",lc=:blue)
        plot!(t,nm(area_type,"WiE"),label="",lc=:red)
        push!(p,p4)
        
        p5 = plot(t,nm(area_type,"VI"),ylabel=latexstring("\\textrm{Mem. Pot.}"),label="",lc=:blue)
        plot!(t,nm(area_type,"VE"),label="",lc=:red)
        push!(p,p5)
        
        responsetype = Bandpass(nm.hp.bandpass[1],nm.hp.bandpass[2],fs=1000/(nm.hp.saveat))
        designmethod = Butterworth(4)
        EEG = filt(digitalfilter(responsetype,designmethod),nm("EEGraw"))
        EEG = reshape(EEG,length(EEG))
        p6 = plot(t,EEG,ylabel=latexstring("\\textrm{EEG}"),lc=:black,label="")
        push!(p,p6)

        per_ = periodogram(EEG,fs = 1000/(nm.hp.saveat))
        p7 = plot(per_.freq,per_.power,label="",ylabel=latexstring("\\textrm{Power}"),xlabel=latexstring("\\omega"),xlims = (0,50))
        push!(p,p7)

        l = @layout [grid(1,1) grid(1,1) grid(1,1);grid(1,1) grid(1,1); grid(1,1) grid(1,1)]
        plot(p...,layout=l;kw...)
        area_name = String(Symbol(typeof(area.pop1).parameters[1]))
        if name !=nothing
            savefig("$(name)-$(area_name).pdf")
        else
            savefig("$(area_name).pdf")
        end
        ctr = ctr + 1
    end
end

function plot_syn(nm::BioNM,fac=1,xlabel="time (ms)", plotname=nothing;kw...)
    name = []
    for area in nm.areas
        type_ = typeof(area.pop1).parameters
        lb = String.(Symbol.(type_))
        if lb[1][1] == 'C' && lb[2][1] == 'E'
            lb = "EEG"
        else
            lb = prod([lb[1][1],",",lb[2][1]])
        end
        push!(name,lb)
        
        type_ = typeof(area.pop2).parameters
        lb = String.(Symbol.(type_))
        if lb[1][1] == 'C' && lb[2][1] == 'E'
            lb = "EEG"
        else
            lb = prod([lb[1][1],",",lb[2][1]])
        end
        push!(name,lb)
    end

    syn_curr = nm("syn_curr")
    ind_ = findall(x->x=="EEG",name)
    
    responsetype = Bandpass(nm.hp.bandpass[1],nm.hp.bandpass[2],fs=1000/(nm.hp.saveat))
    designmethod = Butterworth(4)
    EEG = filt(digitalfilter(responsetype,designmethod),nm("EEGraw"))
    EEG = reshape(EEG,length(EEG))

    syn_curr[:,ind_] = EEG
    p = [] 
    for i in 1:size(syn_curr)[2]
        lab = name[i]
        color = nothing
        if lab[3] == 'E'
            color = :blue
        elseif lab[3] == 'I'
            color = :red
        else
            color = :black
        end

        push!(p,plot(nm.sol.t/fac,syn_curr[:,i],lc=color,label=lab,xlabel=xlabel))
    end

        per_ = periodogram(EEG,fs = 1000/(nm.hp.saveat))
        push!(p,plot(per_.freq,per_.power,label="",ylabel=latexstring("\\textrm{Power}"),xlabel=latexstring("\\omega"),xlims = (0,50)))

    ptitle = plot(title="$(nm.hp.excite[1])pA thalamic stim., between t=$(nm.hp.excite[2]/fac), t=$(nm.hp.excite[3]/fac)",grid=false,axis=nothing,border=:none)
    if name == nothing
        return plot([ptitle; p]...,layout=@layout[A{0.01h}; B C; D E F])
    else
        plot([ptitle; p]...,layout=@layout[A{0.01h}; B C; D E F])
        savefig("$(plotname)-SynapticCurrents.pdf")
    end
end
      
function plot_analysis(nm::BioNM,name::String; kw...)
    responsetype = Bandpass(nm.hp.bandpass[1],nm.hp.bandpass[2],fs=1000/(nm.hp.saveat))
    designmethod = Butterworth(4)
    EEG = filt(digitalfilter(responsetype,designmethod),nm("EEGraw"))
    EEG = reshape(EEG,length(EEG))
    
    p1 = psd(EEG; nfft=3*1000/nm.hp.saveat,xlims=(1,40),window=DSP.hamming(3*1000/nm.hp.saveat),fs=1000/nm.hp.saveat)
    y = tfd(EEG,Spectrogram(nfft=3000/nm.hp.saveat,noverlap=1000/nm.hp.saveat,window=DSP.hamming); fs=1000/nm.hp.saveat)
    p2 = plot(y,clim(0,20),ylims=(0,20))
    plot(p1,p2,layout=(1,2);kw...)
    savefig("$(name)-PSD.pdf")
end

