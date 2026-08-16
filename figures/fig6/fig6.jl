#=

    FIGURE 6
    ========

=#

using RobustAdaptiveTherapy
using Plots, StatsPlots

include("../defaults.jl")

## Setup strategy

Δ = 30.0
α = [0.5,1.0]
tmax = 10000.0
ω = 1e-2
n = 10000

## Patient 99

    met_idx_sample = begin
        Q = sample(patient_fit(99),n)
        sample_metastasis(Q,ω,Δ,α;tmax,expected=false,remove_inf=true)
    end

    met_idx_expected = begin
        Q = sample(patient_fit(99),n)
        sample_metastasis(Q,ω,Δ,α;tmax,expected=true,remove_inf=true)
    end

    fig6a = scatter(reverse(eachcol(met_idx_sample))..., c=:purple, α=0.05, label="")
    fig6b = scatter(reverse(eachcol(met_idx_expected))..., c=:purple, α=0.05, label="")
 
    prop_at_better = count(met_idx_sample[:,2] .> met_idx_sample[:,1]) / n
    prop_at_better = round(prop_at_better, sigdigits=4)
    annotate!(fig6a, 1000, 8000, text("$prop_at_better", :black, :left))

## Joint posterior

    Q = sample_joint_posterior(n)
    met_all_sample = sample_metastasis(Q,ω,Δ,α;tmax,expected=false,remove_inf=false)
    ttp = sample_ttp(Q,Δ,α;tmax,remove_inf=false)

    # Colour based on finite TTP
    idx = isfinite.(ttp[:,1])

    fig6c = scatter(reverse(eachcol(met_all_sample))..., group=idx,c=[:purple :orange], α=0.05, label=["Finite TTP"  "Inf TTP"])

    prop_at_better_all_prog = count(met_all_sample[idx,2] .> met_all_sample[idx,1]) / n
    prop_at_better_all_prog = round(prop_at_better_all_prog, sigdigits=4)
    prop_at_better_all_noprog = count(met_all_sample[.!idx,2] .> met_all_sample[.!idx,1]) / n
    prop_at_better_all_noprog = round(prop_at_better_all_noprog, sigdigits=4)

    annotate!(fig6c, 1000, 8000, text("$prop_at_better_all_prog", :purple, :left))
    annotate!(fig6c, 1000, 7000, text("$prop_at_better_all_noprog", :orange, :left))

## Figure

    fig6 = plot(fig6a,fig6b,fig6c, layout=(1,3), aspect_ratio=:equal,
        xlabel="adapt", ylabel="cont", size=(700,300), 
        xlim=(0.0,10000.0), ylim=(0.0,10000.0))

    [plot!(fig6, subplot=i, x->x, c=:black, ls=:dash, lw=1.5, label="") for i in 1:3]

    savefig(fig6, split("$(@__FILE__).svg",".")[1] * ".svg")
    fig6