#=

    FIGURE 3
    ========

=#

using RobustAdaptiveTherapy
using Plots, StatsPlots
using StatsBase

include("../defaults.jl")


## (a) - Individual estimates for dD
patients = patient_fit_list(model=:strobl)
fig3a = plot()
for i = eachindex(patients)
    res = patient_fit(patients[i];model=:strobl);
    map = get_map(res)[:dD]
    ci95 = quantile(hcat(sample(res,1000)...)[par_idx,:],[0.005,0.995])

    plot!(fig3a, [i,i], ci95, c=:black, label=i == 1 ? "95% CrI" : "")
    scatter!(fig3a, [i], [map], c=:black, label=i == 1 ? "MAP" : "")

end
plot!(fig3a,xticks=(eachindex(patients),patients), xrotation=90,size=(600,200))

## (b) - Joint posterior
X = hcat(sample_joint_posterior(10000)...)

par_idx = findfirst(RobustAdaptiveTherapy.par_names .== "dD")

fig3b = density(X[par_idx,:])

## Figure 3
fig3 = plot(fig3a, fig3b, layout=@layout([a{0.7w} b]), size=(700,200))
savefig(fig3, "$(@__DIR__)/fig3.svg")