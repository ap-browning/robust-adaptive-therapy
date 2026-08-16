#=

    FIGURE 8
    ========

    Stochastic PSA model

    Also produces Fig S4 and S5

=#

using RobustAdaptiveTherapy
using Plots
using Distributions

include("../defaults.jl")

## Patient `idx`

    # Patient 99
    idx = 99

    # Get fits
    res = patient_fit(idx, model=:stochastic_psa);

    # Load data
    day_data, psa_data, Tx_data = patient_data(idx)


## (a,b,c) MAP

    # Get MAP
    q = get_map(idx; model=:stochastic_psa)

    # Simulate trajectories
    psa_sim = [begin
        stoch_psa_simulate_model(q,Tx_data;tmax=maximum(day_data),output=:psa)  
    end for _ = 1:20];

    # Solve for mean and variance
    n, p̂, σ, Σ = stoch_psa_solve_model(q,Tx_data;tmax=maximum(day_data),output=:all);

    # Plot PSA
    fig8a = plot(p̂, ribbon=σ, c=:red, lw=2.0, xlim=(0.0,maximum(day_data)), label="Mean ± Std", fα=0.2)
    plot!(fig8a, psa_sim, label="", c=:black, α=0.2)
    scatter!(fig8a, day_data, psa_data, c=:red, label="Data", ms=4, msw=0.0)
    plot!(fig8a, xlabel="Time [d]", ylabel="Normalised PSA")
    plot!(fig8a, ylim=(-0.2,3.0))

    # Plot n
    plot!(fig8a, t -> n(t) / n(0), c=:blue, lw=2.0, label="n")
    plot!(twinx(), ylim = [ylims(fig8a)...] * n(0), label="", ylabel="Normalised GTV", grid=:off)
    plot!(fig8a, grid=:on)

    # Plot treatment
    fig8b = plot(Tx_data, xlim=extrema(day_data),frange=0.0,lw=0.0,c=:blue,label="",yticks=[],ylabel="Tx")
    plot!(fig8b, ywiden=false)

    # Simulate α
    α_sim = [begin
        sol = stoch_psa_simulate_model(q,Tx_data;tmax=maximum(day_data),output=:raw) 
        t -> sol(t)[4] 
    end for _ = 1:20];

    # Mean and variance
    μ_α = q[1]
    σ_α = q[3] / sqrt(2 * q[2])

    # Plot α
    fig8c = plot(t -> μ_α, ribbon=t -> σ_α, c=:red, label="Mean ± Std", xlim=(0.0,maximum(day_data)), fα=0.3)
    plot!(fig8c, α_sim, label="", c=:black, α=0.2)
    
    plot!(fig8c,xlabel="Time [d]", ylabel="PSA Secretion Rate", xlim=(0.0,maximum(day_data)))
    
    

## Figure 8
fig8 = plot(fig8a, fig8b, fig8c, layout=@layout([a{0.6h}; b{0.1h}; c]), size=(400,400), label="", xwiden=true)
#savefig(fig8, "$(@__DIR__)/fig8.svg")

## Associated supplementary figures

## Fig S5 Posterior Prediction

    T = range(0.0, maximum(day_data), 201)
    p₀_idx = findfirst(res.name_map.parameters .== :p₀)

    # Function to get solution for each sample (PSA)
    func_psa = q -> 
        t -> stoch_psa_solve_model(q,Tx_data;tmax=maximum(day_data),output=:psa)(t)

    # 95% credible interval
    l_psa,u_psa = get_quantiles(func_psa, T, res, 1000)

    # At posterior mode
    m_psa = func_psa(get_map(res))

    figS5 = plot(T, u_psa, frange=l_psa,c=:red, α=0.2, lw=0.0, label="95% CrI (PSA)")
    plot!(figS5, m_psa, c=:red, label="MAP (PSA)", xlim=(0.0,1000.0))

    # GTV (Show with left axis mode unity)
    scale = stoch_psa_solve_model(get_map(res),Tx_data;tmax=maximum(day_data),output=:all)[1](0)
    func_gtv = q -> begin
        sol = stoch_psa_solve_model(q,Tx_data;tmax=maximum(day_data),output=:all)[1]
        t -> sol(t) / scale
    end

    # 95% credible interval
    l_gtv,u_gtv = get_quantiles(func_gtv, T, res, 1000)

    # At posterior mode
    m_gtv = func_gtv(get_map(res))
    
    # Plot
    plot!(figS5, T, u_gtv, frange=l_gtv,c=:blue, α=0.2, lw=0.0, label="95% CrI (GTV)")
    plot!(figS5, m_gtv, c=:blue, label="MAP (GTV)")
    
    # Add data
    scatter!(figS5, day_data, psa_data, c=:red, msw=0.0, label="Data")
    plot!(figS5, ylim=(-0.2,3.0))

    # Add second axis
    plot!(twinx(), ylim = [ylims(figS5)...] * scale, label="", ylabel="GTV (Normalised)", grid=:off)
    plot!(figS5, grid=:on)

    #savefig(figS5, "$(@__DIR__)/figS5.svg")