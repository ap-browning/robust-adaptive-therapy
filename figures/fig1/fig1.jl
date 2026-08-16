#=

    FIGURE 2
    ========

    Demonstrate statistical model (link between PSA and tumour volume, and time-varying parameters)

=#

using Plots
using DifferentialEquations
using Distributions
using RobustAdaptiveTherapy

include("../defaults.jl")

# Use MAP for patient 26
idx = 99
σ,τ = get_map(idx)[[:σ,:τ]]

## (b) Autocorrelation function
acf = δt -> exp(-abs(δt) * log(2) / τ)

fig2b = plot(acf, c=:red, lw=2.0, xlim=(0.0,500.0),
    xlabel="Time [d]", ylabel="PSA autocorrelation", label="", widen=true)

## (c) Sampled PSA trajectory (dual axis)
t = range(0.0,500.0,501)
psa = rand(psa_noise_model(t,σ,τ)) .+ 1.0

fig2c = plot(xlabel="Time [d]", ylabel="Normalised PSA", box=:on)

# Normalised
plot!(fig2c,t,psa / psa[1],c=:red,lw=1.5,label="PSA")
plot!(fig2c,t[1:30:end],(psa / psa[1])[1:30:end],c=:red,m=:circle,msw=1.0,msc=:black,lw=0.0,ms=5,label="Recorded PSA")
hline!([1.0],c=:red,ls=:dash,lw=1.5,label="PSA Baseline")
hline!([1.0 / psa[1]],c=:black,ls=:dash,lw=1.5,label="GTV Baseline")

# Unnormalised
plot!(twinx(), ylim=[ylims(fig2c)...] * psa[1], label="", ylabel="Unnormalised",grid=:off)
plot!(fig2c, grid=:on, box=:on)

## Save
fig2 = plot(fig2b, fig2c, layout=@layout([a b{0.6w}]), size=(600,220), link=:x)
add_plot_labels!(fig2,offset=1)
savefig(fig2, split("$(@__FILE__).svg",".")[1] * ".svg")