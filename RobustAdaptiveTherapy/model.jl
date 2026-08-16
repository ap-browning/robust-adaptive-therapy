#################################
## PARAMETER HANDLING
#################################

function get_pars(q,Tx=default_Tx())
    rS, cR, d̂T, dD, n₀, fR, p₀, σ, τ = q
    rR = rS * cR
    dT = d̂T * rS
    s₀ = (1 - fR) * n₀
    r₀ = fR * n₀
    return (rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx)
end


#################################
## ODE MODEL
#################################

# Strobl model
function f!_strobl(du,u,p,t)
    s,r,drug,∫ = u
    rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = p
    if isa(Tx,Function)
        drug = Tx(t)
    end
    du[1] = rS * s * (1 - (s + r)) * (1 - dD * drug) - dT * s
    du[2] = rR * r * (1 - (s + r)) - dT * r
    du[3] = 0.0
    du[4] = max(0.0,s + r)^(2/3)
end

# Density-independent drug action
function f!_alternative(du,u,p,t)
    s,r,drug = u
    rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = p
    if isa(Tx,Function)
        drug = Tx(t)
    end
    du[1] = rS * s * (1 - (s + r)) * (1 - dD * drug) - dT * s
    #du[2] = rR * r * (1 - (s + r)) * (1 - 0.5 * dD * drug) - dT * r     # 50% less affect of drug than sensitive cells
    du[2] = rR * r * (1 - (s + r)) * (1 - 0.5 * drug) - dT * r     # 50% rediction in rR
    du[3] = 0.0
end

"""
    solve_model(p;tmax=10000.0,output=:psa,model=:strobl)

Solve the model. The tuple p contains both the parameters and information about the treatment schedule. For example:
    `p = get_pars(q)`
which uses (by default) the default treatment strategy.

    solve_model(q,Tx=default_Tx();kwargs...)

Alternative calling sequence.
"""
function solve_model(p::Tuple;tmax=10000.0,output=:psa,model=:strobl)
    # Dispatch to stochastic model (if needed)
    if model == :stochastic
        return stoch_solve_model(p;tmax,output)
    elseif model == :stochastic2
        return stoch2_solve_model(p;tmax,output)
    end

    rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = p
    ode! = eval(Meta.parse("f!_" * String(model)))

    # Solve the ODE
    if isa(Tx,Function)
        # Tx is given as a function of time
        sol = solve(ODEProblem(ode!,[s₀,r₀,1,0],(0.0,tmax),p);verbose=false)
    else
        # Tx is given as a callback (i.e., strategy)
        sol = solve(ODEProblem(ode!,[s₀,r₀,1,0],(0.0,tmax),p);verbose=false,callback=Tx)
    end
    psa = t -> sum(sol(t)[1:2]) / (s₀ + r₀)
    drug = t -> sol(t)[3] > 0.5
    rprp = t -> sol(t)[2] / sum(sol(t)[1:2])
    if output == :raw
        return sol
    elseif output == :psa
        return psa
    else
        return psa, drug, rprp
    end
end
function solve_model(q,Tx=default_Tx();model=:strobl,kwargs...)
    # Dispatch to stochastic model (if needed)
    if model == :stochastic
        return stoch_solve_model(stoch_get_pars(q,Tx);kwargs...)
    elseif model == :stochastic2
        return stoch2_solve_model(stoch_get_pars(q,Tx);kwargs...)
    end
    return solve_model(get_pars(q,Tx);kwargs...)
end


"""
    create_Tx(α,Δ,ε=0.0;tmax=10000.0)
    
Create a treatment strategy callback, based on a decision interval of Δ days. Drug is turns off if psa(t) drops below α[1] * psa(0.0), on if rises above α[2] * psa(0.0). Optionally include a vector `ε` representing (potentially correlated) measurement errors. 
"""
function create_Tx(α,Δ,ε=0.0;tmax=10000.0,model=:strobl)
    return PresetTimeCallback(0:Δ:tmax, int -> begin
        if model == :stochastic
            s,r,K,drug = u
            KT,Kξ,Kη,rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = int.p
        elseif model == :stochastic_psa
            s,r,p,α,drug = u
            α₀,α̂,ξα,ηα,β,rS,rR,dT,dD,s₀,r₀,p₀,σ,Tx = int.p
        else
            rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = int.p
            s,r,drug = int.u
        end  

        # Calculate PSA mean (i.e., normalised cell count)
        psa = sum(s + r) / (s₀ + r₀)
        # Add noise, if required, renormalise
        if isa(ε,Vector)
            psa = (psa + ε[Int(int.t / Δ + 1)]) / (1 + ε[1])    
        end
        # Make a decision about the drug
        if psa < α[1]
            drug = 0
        elseif psa ≥ α[2]
            drug = 1
        end
        if model != :stochastic
            int.u[3] = drug
            
        else
            int.u[4] = drug
        end
    end)
end

"""
    The default strategy is that of Strobl. et al, with a daily decision interval.
"""
default_Tx(ε=0.0;kwargs...) = create_Tx([0.5,1.0],7.0,ε;kwargs...)

"""
    Correlated PSA noise model
"""
function psa_noise_model(t,σ::Number,τ::Number)
    acf = δt -> σ^2 * exp(-abs(δt) * log(2) / τ)
    MvNormal([acf(t₁ - t₂) for t₁ in t, t₂ in t])
end
psa_noise_model(t,q) = psa_noise_model(t,q[end-1:end]...)


#################################
## DIAGNOSTICS, ANALYSIS
#################################

function model_all_equilibria(q;drug=true,stability=true)
    rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = get_pars(q)
    if drug
        rS = rS * (1 - dD)
    end
    x₁ = [0,0]
    x₂ = [(rS - dT) / rR,0]
    x₃ = [0,(rR - dT) / rR]
    if !stability
        return [x₁,x₂,x₃]
    end
    function isstable(x)
        s,r = x
        jac = [-dT + rS * (1 - r - s) - rS * s -rS * s;
               -r * rR  -dT - r * rR + rR * (1 - r - s)]
        all(real.(eigvals(jac)) .< 0.0) 
    end
    # s₁ = all([rR - dT; rS - dT] .< 0.0)
    # s₂ = all([dT - rS; rR - rS] .< 0.0)
    # s₃ = all([dT - rR; rS - rR] .< 0.0)
    s₁ = isstable(x₁)
    s₂ = isstable(x₂)
    s₃ = isstable(x₃ )
    return ([x₁,x₂,x₃],[s₁,s₂,s₃])
end
function model_stable_equilibria(q;drug=true)
    X,S = model_all_equilibria(q;drug,stability=true)
    X[S][1]
end
function model_patient_map_progress(idx)
    q = get_map(idx)
    model_stable_equilibria(q)[2] > 0.0
end

function find_ttp(psa,drug,β=1.2;tmax=10000.0)
    T = 0:1.0:tmax
    I = (psa.(T) .> β) .& drug.(T)
    if any(I)
        idx = findlast(I .== false) + 1
        if idx == length(T) + 1
            return Inf
        end
        t₂ = T[findlast(I .== false) + 1]
        # Setup bracketing interval...
        t₁ = max(0.0,t₂ - 1.0)
        func = t -> psa(t) - β
        if sign(t₁) != sign(t₂)
            return find_zero(func, [t₁,t₂])
        else
            t₁ = T[findlast(psa.(T) .< β)]
            return find_zero(func, [t₁,t₂])
        end
    else
        return Inf
    end
end

function sample_ttp(Q,Δ,α;tmax=10000.0,model=:strobl,remove_inf=true,kwargs...)

    # Treatment strategies
    Tx_int = ε -> create_Tx(α,Δ,isa(ε,Distribution) ? rand(ε) : ε;tmax,model)
    Tx_con = t -> true

    # PSA noise model
    ε = q -> psa_noise_model(0:Δ:tmax,q)

    function get_psa_drug(q,ε=nothing)
        if isnothing(ε)
            return solve_model(q,Tx_con;tmax,output=:all,model,kwargs...)[1:2]
        else
            return solve_model(q,Tx_int(ε);tmax,output=:all,model,kwargs...)
        end
    end
    function get_ttp(q,ε=nothing)
        psa,drug = get_psa_drug(q,ε)
        find_ttp(psa,drug;tmax)
    end

    # Store results in a matrix...
    ttp = zeros(length(Q),2)
    @threads for i = eachindex(Q)
        ttp[i,1] = get_ttp(Q[i])
        if isinf(ttp[i,1])
            ttp[i,2] = Inf
        else
            ttp[i,2] = get_ttp(Q[i],ε(Q[i])) 
        end
    end 

    if remove_inf
        ttp = ttp[all(isfinite.(ttp),dims=2)[:],:]
    end

    ttp

end
sample_ttp_joint_posterior(n::Number,args...;kwargs...) = 
    sample_ttp(sample_joint_posterior(n),args...;kwargs...)

function sample_metastasis(Q,ω,Δ,α;tmax=10000.0,model=:strobl,expected=false,remove_inf=false,dt=1.0,kwargs...)

    # Treatment strategies
    Tx_int = ε -> create_Tx(α,Δ,isa(ε,Distribution) ? rand(ε) : ε;tmax,model)
    Tx_con = t -> true

    # PSA noise model
    ε = q -> psa_noise_model(0:Δ:tmax,q)

    function get_sol(q,ε=nothing)
        if isnothing(ε)
            sol = solve_model(q,Tx_con;tmax,output=:raw,model,kwargs...)
        else
            sol = solve_model(q,Tx_int(ε);tmax,output=:raw,model,kwargs...)
        end

    end

    function get_metastasis_time(q,ε=nothing)
        sol = get_sol(q,ε)
        if expected
            ∫n = t -> sol(t)[4]
            func = t -> exp(-ω * ∫n(t))
            if func(10000) > 1e-3
                @warn "Expected metastasis time may be inaccurate; consider increasing tmax"
            end
            return quadgk(func,0.0,tmax)[1]
        else
            n = t -> sum(sol(t)[1:2])
            t = 0.0:dt:tmax
            e = rand.(Poisson.(max.(0.0,ω * n.(t) * dt))) .> 0.0
            return any(e) ? t[findfirst(e .== true)] + rand() * dt : Inf
        end
    end

    # Store results in a matrix...
    meta = zeros(length(Q),2)
    @threads for i = eachindex(Q)
        meta[i,1] = get_metastasis_time(Q[i])
        meta[i,2] = get_metastasis_time(Q[i],ε(Q[i])) 
    end 

    if remove_inf
        meta = meta[all(isfinite.(meta),dims=2)[:],:]
    end

    meta

end
sample_metastasis_joint_posterior(n::Number,args...;kwargs...) = 
    sample_metastasis(sample_joint_posterior(n),args...;kwargs...)