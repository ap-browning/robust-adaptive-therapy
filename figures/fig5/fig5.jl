#=

    FIGURE 5
    ========

=#

using RobustAdaptiveTherapy
using Plots, StatsPlots
using JLD2

include("../defaults.jl")

## Setup various strategies
Δ = [7,30,60]
α = [[0.2,0.8],[0.2,1.0],[0.5,0.8],[0.5,1.0]]
tmax = 10000.0
n = 10000

@time ttp = [sample_ttp_joint_posterior(n,Δᵢ,αᵢ;tmax,remove_inf=false) for αᵢ = α, Δᵢ = Δ]
@save "$(@__DIR__)/ttp.jld2" ttp
@load "$(@__DIR__)/ttp.jld2" ttp

prop_noprogress = count(isinf.(vcat(ttp[:]...)[:,1])) / (n * prod(size(ttp)))

## (a) Scatter plot for distribution in figure 3

ttp_fig3 = ttp[findfirst(α .== [[0.5,1.0]]),findfirst(Δ .== 30)]
ttp_fig3_finite = ttp_fig3[isfinite.(ttp_fig3[:,1] .* ttp_fig3[:,2]),:]

fig5a = scatter(reverse(eachcol(ttp_fig3_finite))..., c=:purple, α=0.2, label="")
plot!(fig5a, x->x, c=:black, ls=:dash, lw=1.5, label="",
        aspect_ratio=:equal,xlim=(0.0,3000.0),ylim=(0.0,3000.0))
plot!(fig5a, xlabel="TTP (adaptive)", ylabel="TTP (continuous)")

# Annotations
prop_progress = count(isfinite.(vcat(ttp...)[:,1])) / length(vcat(ttp...)[:,1])
prop_progress = round(prop_progress, sigdigits=4)
prop_improve = count(ttp_fig3_finite[:,1] .< ttp_fig3_finite[:,2]) / size(ttp_fig3_finite,1)
prop_improve = round(prop_improve, sigdigits=4)
annotate!(fig5a,(50.0,2800.0,("p (progress) = $prop_progress \np (improve) = $prop_improve",:left,10)))

## (b) Violin plots of the improvement

improve = [filter(isfinite, ttpᵢ[:,2] .- ttpᵢ[:,1]) for ttpᵢ in ttp]
h = 0.3
#clip = 0.05
c = palette(:Set2_3)

fig5b = plot()
for i = axes(improve)[1], j = axes(improve)[2]
    y = improve[i,j]
    boxplot!(fig5b,fill(i + (-h:h:h)[j],size(y)), y, bar_width=0.1, c=c[j], 
        label=i == 1 ? "Δ = $(Δ[j]) d" : "", lw=1.0)
end
plot!(fig5b, ylim=[-500.0,1500.0], widen=true, xticks=(eachindex(α), α), 
    ylabel="Adaptive Improvement", xlabel="α")

## Figure 5
fig5 = plot(fig5a, fig5b, size=(900,300), layout=@layout([a b{0.6w}]))
savefig(fig5, split("$(@__FILE__).svg",".")[1] * ".svg")
fig5