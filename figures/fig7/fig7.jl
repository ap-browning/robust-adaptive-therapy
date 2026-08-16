#=

    FIGURE 7
    ========

    Stochastic models

=#

using RobustAdaptiveTherapy
using Plots
using Distributions

include("../defaults.jl")

## Load data

    # Patient 2
    idx = 2

    # Get fits
    res = patient_fit(idx, model=:stochastic)

    # Load data
    day_data, psa_data, Tx_data = patient_data(idx)


## Simulate

    # Get MAP
    q = sample(res,1)[1]

    # Simulate trajectories
    sde_sim = [begin
        stoch_solve_model(q,Tx_data;tmax=maximum(day_data),output=:psa)  
    end for _ = 1:20];

    # Plot PSA
    fig7a = plot(sde_sim, xlim=extrema(day_data), label="", c=:black, α=0.2)
    scatter!(fig7a, day_data, psa_data, c=:red, label="Data", ms=4, msw=0.0)

    # Plot treatment
    fig7b = plot(Tx_data, xlim=extrema(day_data),frange=0.0,lw=0.0,c=:blue,label="",yticks=[],ylabel="Tx")
    plot!(fig7b, ywiden=false)

    # Simulate K
    K_sim = [begin
        sol = stoch_solve_model(q,Tx_data;tmax=maximum(day_data),output=:raw) 
        t -> sol(t)[3]
    end for _ = 1:20];

    # Plot K
    fig7c = plot(K_sim, xlim=(0.0,maximum(day_data)), label="", c=:black, α=0.2)
    plot!(fig7c,xlabel="Time [d]", ylabel="K")

# Figure 7
fig7 = plot(fig7a, fig7b, fig7c, layout=@layout([a{0.6h}; b{0.1h}; c]), size=(400,400), label="", xwiden=true)
savefig(fig7, "$(@__DIR__)/fig7.svg")
fig7
