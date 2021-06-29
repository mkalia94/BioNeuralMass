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

    p1 = plot(t,nm("NaCiI"),xlabel=latexstring("[\\textrm{Na}^+]"),lc=:blue)
    plot!(t,nm("NaCiE"),lc=:red)
    plot!(t,nm("NaCe"),lc=:green)
    push!(p,p1)

    p2 = plot(t,nm("KCiI"),xlabel=latexstring("[\\textrm{K}^+]"),lc=:blue)
    plot!(t,nm("KCiE"),lc=:red)
    plot!(t,nm("KCe"),lc=:green)
    push!(p,p2)

    p3 = plot(t,nm("ClCiI"),xlabel=latexstring("[\\textrm{Cl}^-]"),lc=:blue)
    plot!(t,nm("ClCiE"),lc=:red)
    plot!(t,nm("ClCe"),lc=:green)
    push!(p,p3)

    p4 = plot(t,nm("WiI"),xlabel=latexstring("\\textrm{Vol.}"),lc=:blue)
    plot!(t,nm("WiE"),lc=:red)
    plot!(t,nm("We"),lc=:green)
    push!(p,p4)
    
    p5 = plot(t,nm("V_I"),xlabel=latexstring("\\textrm{Mem. Pot.}"),lc=:blue)
    plot!(t,nm("V_E"),lc=:red)
    push!(p,p5)
    
    responsetype = Bandpass(nm.hp.bandpass[1],nm.hp.bandpass[2],fs=1/(nm.hp.saveat*60))
    designmethod = Butterworth(4)
    EEG = filt(digitalfilter(responsetype,designmethod),nm("I_E"))
    p6 = plot(t,EEG,xlabel=latexstring("\\textrm{EEG}"),lc=:black)
    push!(p,p6)
    
    l = @layout [grid(1,1) grid(1,1) grid(1,1);grid(1,1) grid(1,1); grid(1,1)]
    plot(p...,layout=l;kw...)
end

