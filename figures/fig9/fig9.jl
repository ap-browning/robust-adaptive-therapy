#=

    FIGURE 9
    ========

=#

using RobustAdaptiveTherapy
using Plots
using Distributions
using JLD2
using Random

include("../defaults.jl")

## ALTERNATIVE MODEL PARAMETERS

    # Get random parameter from joint posterior
    Random.seed!(4)
    q = sample_joint_posterior(model=:strobl)[1]
    Random.seed!()

# SETUP TREATMENT SCHEDULE

    tmax = 500.0       # Maximum time (such that TTP is always less than tmax)
    Δ = 14              # Decision interval
    day = 0:Δ:tmax      # Observation times (for decisions)
    α = [0.5,1.0]       # Off when Rel. PSA drops below α[1], on when above α[2]

    # Adaptive therapy
    Tx_at = ε -> create_Tx(α,Δ,isa(ε,Distribution) ? rand(ε) : ε;tmax)

    # Continuous therapy
    Tx_con = t -> true 

    # Intermittent therapy
    Tx_int = t -> mod(t, 200) < 100     # Intermittent

    # ε for adaptive therapy (data)
    ε_sample = rand(psa_noise_model(day,q))

## SYNTHETIC DATA FROM ALTERNATIVE MODEL

    # Solve model
    sol_psa, Tx_data, rprp = solve_model(q, Tx_at(ε_sample);tmax, output=:all, model=:alternative);

    # Plot measured PSA
    psa_data = sol_psa.(day) + ε_sample
    psa_data = psa_data / psa_data[1]

    p1 = plot(sol_psa, xlim=(0.0,tmax))
    scatter!(day, psa_data)


## FIT TO STROBL MODEL

    #res = patient_fit(day, psa_data, Tx_data; model=:strobl, tryagain=20);
    #@save "$(@__DIR__)/fig9.jld2" res q ε_sample psa_data
    @load "$(@__DIR__)/fig9.jld2" res q ε_sample psa_data

## PLOTS

    T = range(0.0,tmax, 201)
    Tex = range(0.0, 4*tmax)    # For extrapolation

    n = 1000 # Number of samples for credible intervals

    p₀_idx = findfirst(res.name_map.parameters .== :p₀)

    # (a) - Fit
    
        # Function of parameters to plot
        func = q -> t -> q[p₀_idx] * solve_model(q, Tx_data; tmax, output=:psa, model=:strobl)(t)
    
        # 95% credible interval
        l,u = get_quantiles(func, T, res, n)

        # At posterior mode
        m = func(get_map(res))

        fig9a = plot(T, u, frange=l,c=:blue, α=0.2, lw=0.0, label="95% CrI")
        plot!(fig9a, m, c=:blue, label="MAP (Misspecified)")
        plot!(fig9a, sol_psa, c=:red, label="True")
        scatter!(fig9a, day,psa_data, c=:red, msw=0.0, label="Data")
        plot!(fig9a, ylim=(0.0,:auto), widen=true)

        ###########
        ## Plot resistant population
        sol_true = t -> solve_model(q, Tx_data; tmax=maximum(Tex), output=:raw, model=:alternative)(t)[2];

        # Function of parameters to plot
        func = q -> t -> solve_model(q, Tx_data; tmax=maximum(Tex), output=:raw, model=:strobl)(t)[2]

        # 95% credible interval
        l,u = get_quantiles(func, T, res, n)

        # At posterior mode
        m = func(get_map(res))

        fig9d = plot(T, u, frange=l,c=:blue, α=0.2, lw=0.0, label="95% CrI")
        plot!(fig9d, m, c=:blue, label="MAP (Misspecified)")
        plot!(fig9d, sol_true, c=:red, label="True")

    # (b) - Continuous
    
        # Truth
        sol_psa_cont_true = solve_model(q, Tx_con; tmax=maximum(Tex), output=:psa, model=:alternative)

        # Function of parameters to plot
        func = q -> t -> solve_model(q, Tx_con; tmax=maximum(Tex), output=:psa, model=:strobl)(t)

        # 95% credible interval
        l,u = get_quantiles(func, Tex, res, n)

        # At posterior mode
        m = func(get_map(res))

        fig9b = plot(Tex, u, frange=l,c=:blue, α=0.2, lw=0.0, label="95% CrI")
        plot!(fig9b, m, c=:blue, label="MAP (Misspecified)")
        plot!(fig9b, sol_psa_cont_true, c=:red, label="True")

        ###########
        ## Plot resistant population
        sol_true = t -> solve_model(q, Tx_con; tmax=maximum(Tex), output=:raw, model=:alternative)(t)[2];

        # Function of parameters to plot
        func = q -> t -> solve_model(q, Tx_con; tmax=maximum(Tex), output=:raw, model=:strobl)(t)[2]

        # 95% credible interval
        l,u = get_quantiles(func, Tex, res, n)

        # At posterior mode
        m = func(get_map(res))

        fig9e = plot(Tex, u, frange=l,c=:blue, α=0.2, lw=0.0, label="95% CrI")
        plot!(fig9e, m, c=:blue, label="MAP (Misspecified)")
        plot!(fig9e, sol_true, c=:red, label="True")
        plot!(fig9e, ylim=(0.0,1.0), widen=true)

    # (c) - Intermittent
    
        # Truth
        sol_psa_int_true = solve_model(q, Tx_int; tmax=maximum(Tex), output=:psa, model=:alternative)

        # Function of parameters to plot
        func = q -> t -> solve_model(q, Tx_int; tmax=maximum(Tex), output=:psa, model=:strobl)(t)

        # 95% credible interval
        l,u = get_quantiles(func, Tex, res, n)

        # At posterior mode
        m = func(get_map(res))

        fig9c = plot(Tex, u, frange=l,c=:blue, α=0.2, lw=0.0, label="95% CrI")
        plot!(fig9c, m, c=:blue, label="MAP (Misspecified)")
        plot!(fig9c, sol_psa_int_true, c=:red, label="True")

        ###########
        ## Plot resistant population
        sol_true = t -> solve_model(q, Tx_int; tmax=maximum(Tex), output=:raw, model=:alternative)(t)[2];

        # Function of parameters to plot
        func = q -> t -> solve_model(q, Tx_int; tmax=maximum(Tex), output=:raw, model=:strobl)(t)[2]

        # 95% credible interval
        l,u = get_quantiles(func, Tex, res, n)

        # At posterior mode
        m = func(get_map(res))

        fig9f = plot(Tex, u, frange=l,c=:blue, α=0.2, lw=0.0, label="95% CrI")
        plot!(fig9f, m, c=:blue, label="MAP (Misspecified)")
        plot!(fig9f, sol_true, c=:red, label="True")
        plot!(fig9f, ylim=(0.0,1.0), widen=true)

# FIGURE 9

fig9 = plot(fig9a, fig9b, fig9c, layout=grid(1,3), size=(800,200), xlabel="Time [d]", ylabel="Normalised PSA")
add_plot_labels!(fig9)

savefig(fig9, "$(@__DIR__)/fig9.svg")
fig9