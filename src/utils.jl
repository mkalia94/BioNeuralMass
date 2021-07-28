import DifferentialEquations: solve
import Plots: plot

function solve(nm::BioNM,hp::HyperParam,areas::Vector{NeuralArea},alg=Tsit5();kw...)
    X0 = []
    for area in areas
        push!(X0,area.X0)
    end
    hp.X0 = vcat(X0...)
    prob = ODEProblem((dx,x,p,t)->nm(areas,hp,dx,x,p,t),hp.X0,[0.0f0,hp.tfinal]) 
    sol = solve(prob,alg;kw...)
    nm.sol = sol
end

function plot(hp::HyperParam,nm::BioNM,areas::Vector{NeuralArea},name=nothing;kw...)
    ctr = 1
    for area in areas
        t = nm.sol.t
        p = []

        p1 = plot(t,nm(areas,hp,ctr,"NaCiI"),ylabel=latexstring("[\\textrm{Na}^+]"),label="",lc=:blue)
        plot!(t,nm(areas,hp,ctr,"NaCiE"),label="",lc=:red)
        push!(p,p1)

        p2 = plot(t,nm(areas,hp,ctr,"KCiI"),ylabel=latexstring("[\\textrm{K}^+]"),label="",lc=:blue)
        plot!(t,nm(areas,hp,ctr,"KCiE"),label="",lc=:red)
        push!(p,p2)

        p3 = plot(t,nm(areas,hp,ctr,"ClCiI"),ylabel=latexstring("[\\textrm{Cl}^-]"),label="",lc=:blue)
        plot!(t,nm(areas,hp,ctr,"ClCiE"),label="",lc=:red)
        push!(p,p3)

        p4 = plot(t,nm(areas,hp,ctr,"WiI"),ylabel=latexstring("\\textrm{Vol.}"),label="",lc=:blue)
        plot!(t,nm(areas,hp,ctr,"WiE"),label="",lc=:red)
        push!(p,p4)
        
        p5 = plot(t,nm(areas,hp,ctr,"VI"),ylabel=latexstring("\\textrm{Mem. Pot.}"),label="",lc=:blue)
        plot!(t,nm(areas,hp,ctr,"VE"),label="",lc=:red)
        push!(p,p5)
        
        responsetype = Bandpass(hp.bandpass[1],hp.bandpass[2],fs=1000/(hp.saveat))
        designmethod = Butterworth(4)
        EEG = filt(digitalfilter(responsetype,designmethod),nm(areas,hp,"EEGraw"))
        EEG = reshape(EEG,length(EEG))
        p6 = plot(t,EEG,ylabel=latexstring("\\textrm{EEG}"),lc=:black,label="")
        push!(p,p6)

        per_ = periodogram(EEG,fs = 1000/(hp.saveat))
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

