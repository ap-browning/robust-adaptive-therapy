#=

    FIGURE 4
    ========

    Looks at TTP distribution variability for patient idx = 99

    Decision interval: Δ = 30 d
    adaptive: α = [0.2,1.0]

=#

using RobustAdaptiveTherapy
using Distributions
using Plots, StatsPlots

include("../defaults.jl")

# Patient `idx` in (a) and (b)
idx = 99
res = patient_fit(idx)

## SETUP TREATMENT SCHEDULE

    tmax = 10000.0      # Maximum time (such that TTP is always less than tmax)
    Δ = 30              # Decision interval
    tobs = 0:Δ:tmax     # Observation times (for decisions)
    α = [0.5,1.0]       # Off when Rel. PSA drops below α[1], on when above α[2]

    Tv_at = ε -> create_Tx(α,Δ,isa(ε,Distribution) ? rand(ε) : ε;tmax)
    Tx_con = t -> true

    ε = q -> psa_noise_model(0:Δ:tmax,q)

    function get_psa_drug(q,ε=nothing)
        if isnothing(ε)
            return solve_model(q,Tx_con;tmax,output=:all)[1:2]
        else
            return solve_model(q,Tv_at(ε);tmax,output=:all)
        end
    end
    function get_ttp(q,ε=nothing)
        psa,drug = get_psa_drug(q,ε)
        find_ttp(psa,drug;tmax)
    end

## (a) - PSA trajectories

    panel1 = begin

        q = get_map(idx)

        # Get continuous therapy TTP
        ε_conᵢ = rand(ε(q))
        psa_con,drug_con = get_psa_drug(q);
        ttp_con = get_ttp(q)

        # Ensemble plot
        ε_intᵢ = rand(ε(q))
        psa_adapt,drug_adapt = get_psa_drug(q,ε_intᵢ)
        ttp_adapt = get_ttp(q,ε_intᵢ)

        # Plots...
        plt1 = hline(α,lw=1.5,c=:red,ls=:dash,label="Target")
        plt2 = hline(α,lw=1.5,c=:red,ls=:dash,label="Target")
        plot!(plt1,psa_con,xlim=(0.0,tmax),c=:blue,lw=1.5,label="Cont.")
        plot!(plt1,psa_adapt,c=:red,lw=1.5,label="Adapt.")
        
        plot!(plt2,tobs,max.(0.0,(psa_con.(tobs) + ε_conᵢ) / (1 + ε_conᵢ[1])),c=:blue,m=:circle,msw=0.0,lw=1.5,label="Cont.")
        plot!(plt2,tobs,max.(0.0,(psa_adapt.(tobs) + ε_intᵢ) / (1 + ε_intᵢ[1])),c=:red,m=:circle,msw=0.0,lw=1.5,label="Adapt.")

        [hline!(plt,[1.2],c=:black,lw=1.5,ls=:dash,label="Progression") for plt = [plt1,plt2]]
    
        plot(plt1,plt2,layout=grid(2,1),ylim=(0.0,3.2),xlim=(0.0,2500.0))

    end

    plot!(panel1,xlabel="Time [d]",ylabel="Norm. PSA",widen=true)

## (b) - TTP distribution

    # Number of simulations
    n = 500

    ## Fixed parameters
    panel2a = begin

        # Continuous strategy
        ttp_con = get_ttp(q)

        # Exact measurements, adaptive strategy
        ttp_adapt_exact = get_ttp(q,0.0*rand(ε(q)))

        # Distribution of TTP
        ttp1 = [get_ttp(q,rand(ε(q))) for _ = 1:n]

        # Plot
        plt = histogram(ttp1,c=:red,lw=0.0,label="Adapt.",ywiden=false)
        vline!([ttp_con],c=:blue,lw=1.5,label="Cont.")
        vline!([ttp_adapt_exact],c=:black,ls=:dash,lw=1.5,label="adapt. (Exact PSA)")

    end

    # Parameters from posterior
    panel2b,panel2c = begin

        # Continuous
        Q = sample(res,n)

        # Distribution of TTP (continuous)
        ttp2_con = [get_ttp(q) for q = Q]

        # Distribution of TTP (adaptive)
        ttp2_int = [get_ttp(q,rand(ε(q))) for q = Q]

        # Plot
        plt1 = plot()
        histogram!(plt1,ttp2_con,c=:blue,lw=0.0,label="Cont.",ywiden=false, ylim=(0.0,89.0))
        histogram!(plt1,ttp2_int,c=:red,lw=0.0,label="Adapt.", binwidth=20, ylim=(0.0,120.0))

        plt2 = scatter(ttp2_int,ttp2_con,c=:purple,α=0.5,label="",
            xlabel="TTP (adaptive)", ylabel="TTP (continuous)")

        (plt1,plt2)

    end

    panel2 = plot(panel2a,panel2b,layout=grid(2,1), link=:x, 
        xlim=(1000.0,2500.0),xlabel="TTP [d]")
    panel3 = plot(panel2c)
    plot!(panel3, x->x, c=:black, ls=:dash, lw=1.5, label="",
        aspect_ratio=:equal,xlim=(1000.0,3000.0),ylim=(1000.0,3000.0))

## Figure 3

fig4 = plot(panel1, panel2, plot!(panel3, aspect_ratio=:equal), size=(900,300),
    layout=@layout([a{0.4w} b c{0.37w}]), link=:x)

savefig(fig4, split("$(@__FILE__).svg",".")[1] * ".svg")
fig4