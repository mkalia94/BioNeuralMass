import DifferentialEquations: solve
import Plots: plot

function solve(nm::BioNM,alg=Tsit5();kw...)
    prob = ODEProblem((x,p,t)->nm(x,p,t),nm.X0,[0.0f0,nm.hp.tfinal]) 
    sol = solve(prob,alg;kw...)
    nm.sol = sol
end

function plot(nm::BioNM;kw...)
    t = nm.sol.t
    p = []

    p1 = plot(t,nm("NaCiI"),ylabel=latexstring("[\\textrm{Na}^+]"),label="",lc=:blue)
    plot!(t,nm("NaCiE"),label="",lc=:red)
    plot!(t,nm("NaCe"),label="",lc=:green)
    push!(p,p1)

    p2 = plot(t,nm("KCiI"),ylabel=latexstring("[\\textrm{K}^+]"),label="",lc=:blue)
    plot!(t,nm("KCiE"),label="",lc=:red)
    plot!(t,nm("KCe"),label="",lc=:green)
    push!(p,p2)

    p3 = plot(t,nm("ClCiI"),ylabel=latexstring("[\\textrm{Cl}^-]"),label="",lc=:blue)
    plot!(t,nm("ClCiE"),label="",lc=:red)
    plot!(t,nm("ClCe"),label="",lc=:green)
    push!(p,p3)

    p4 = plot(t,nm("WiI"),ylabel=latexstring("\\textrm{Vol.}"),label="",lc=:blue)
    plot!(t,nm("WiE"),label="",lc=:red)
    plot!(t,nm("We"),label="",lc=:green)
    push!(p,p4)
    
    p5 = plot(t,nm("V_I"),ylabel=latexstring("\\textrm{Mem. Pot.}"),label="",lc=:blue)
    plot!(t,nm("V_E"),label="",lc=:red)
    push!(p,p5)
    
    responsetype = Bandpass(nm.hp.bandpass[1],nm.hp.bandpass[2],fs=1000/(nm.hp.saveat))
    designmethod = Butterworth(4)
    EEG = filt(digitalfilter(responsetype,designmethod),nm("I_E"))
    p6 = plot(t,EEG,ylabel=latexstring("\\textrm{EEG}"),lc=:black,label="")
    push!(p,p6)

    per_ = periodogram(EEG,fs = 1000/(nm.hp.saveat))
    p7 = plot(per_.freq,per_.power,label="",ylabel=latexstring("\\textrm{Power}"),xlabel=latexstring("\\omega"),xlims = (0,80))
    push!(p,p7)

    l = @layout [grid(1,1) grid(1,1) grid(1,1);grid(1,1) grid(1,1); grid(1,1) grid(1,1)]
    plot(p...,layout=l;kw...)
end

