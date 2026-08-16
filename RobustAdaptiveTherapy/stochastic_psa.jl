#################################
## PARAMETER HANDLING
#################################

"""
   `get_pars`` for the stochastic model
"""
function stoch_psa_get_pars(q,Tx=default_Tx())
    α̂,ξα,ηα,β,rS, cR, d̂T, dD, n₀, fR, p₀, σ = q
    rR = rS * cR
    dT = d̂T * rS
    s₀ = (1 - fR) * n₀
    r₀ = fR * n₀
    return (α̂,ξα,ηα,β,rS,rR,dT,dD,s₀,r₀,p₀,σ,Tx)
end

"""
    `prior` for the stochastic model
"""
prior_α̂ = Uniform(0.0,0.2)
prior_ξα = log(2) * InverseGamma(5,1/100)
prior_ηα = LogUniform(1e-4,1e-2)
prior_β = Uniform(0.0,0.1)
stoch_psa_prior = Product([[
    prior_α̂, prior_ξα, prior_ηα, prior_β]; 
    prior.v[1:end-1]])

stoch_psa_par_names = ["α̂";"ξα";"ηα";"β";par_names[1:end-1]]

"""
    Convert from deterministic to stochastic parameters
"""
function stoch_psa_convert_qode(qode)
    [mean(stoch_psa_prior)[1:4]; qode[1:end-1]]
end

#################################
## SDE MODEL
#################################

function f!_stoch_psa(du,u,par,t)
    s,r,p,α,drug = u
    n = s + r
    α̂,ξα,ηα,β,rS,rR,dT,dD,s₀,r₀,p₀,σ,Tx = par
    if isa(Tx,Function)
        drug = Tx(t)
    end
    du[1] = rS * s * (1 - n) * (1 - dD * drug) - dT * s
    du[2] = rR * r * (1 - n) - dT * r
    du[3] = α * n - β * p
    du[4] = -ξα * (α - α̂)
    du[5] = 0.0
end
function g!_stoch_psa(du,u,par,t)
    du .= 0.0
    α̂,ξα,ηα,β,rS,rR,dT,dD,s₀,r₀,p₀,σ,Tx = par
    du[4] = ηα
end


"""
    stoch_psa_solve_model(p;tmax=10000.0,output=:psa,model=:strobl)

Solve the model. The tuple p contains both the parameters and information about the treatment schedule. For example:
    `p = get_pars(q)`
which uses (by default) the default treatment strategy.


    solve_model(q,Tx=default_Tx();kwargs...)

Alternative calling sequence.
"""
function stoch_psa_simulate_model(p::Tuple;tmax=10000.0,dt=0.5,output=:psa)
    α̂,ξα,ηα,β,rS,rR,dT,dD,s₀,r₀,p₀,σ,Tx = p
    # Sample initial PSA secration rate from stationary distribution 
    α₀ = rand(Normal(α̂, ηα / sqrt(2ξα)))
    # Solve the ODE
    if isa(Tx,Function)
        # Tx is given as a function of time
        sol = solve(SDEProblem(f!_stoch_psa,g!_stoch_psa,[s₀,r₀,p₀,α₀,1],(0.0,tmax),p),EM(),dt=dt;verbose=false)
    else
        # Tx is given as a callback (i.e., strategy)
        sol = solve(SDEProblem(f!_stoch_psa,g!_stoch_psa,[s₀,r₀,p₀,α₀,1],(0.0,tmax),p),EM(),dt=dt;verbose=false,callback=Tx)
    end
    psa = t -> sol(t)[3] / p₀
    drug = t -> sol(t)[5] .> 0.5
    rprp = t -> sol(t)[2] / sum(sol(t)[1:2])
    if output == :raw
        return sol
    elseif output == :psa
        return psa
    else
        return psa, drug, rprp
    end
end
stoch_psa_simulate_model(q,Tx=default_Tx();kwargs...) = stoch_psa_simulate_model(stoch_psa_get_pars(q,Tx);kwargs...)

#################################
## ODE MODEL
#################################

function f!_stoch_psa_ode(du,u,par,t)
    s,r,p̂,A,B,C,D,drug = u
    n = s + r
    α̂,ξα,ηα,β,rS,rR,dT,dD,s₀,r₀,p₀,σ,Tx = par
    if isa(Tx,Function)
        drug = Tx(t)
    else
        error("Derivation of stochastic PSA ODE model does not apply if Tx is a callback!")
    end
    du[1] = rS * s * (1 - n) * (1 - dD * drug) - dT * s
    du[2] = rR * r * (1 - n) - dT * r
    du[3] = α̂ * n - β * p̂
    du[4] = n * exp(β * t) * exp(ξα * t)
    du[5] = n * exp(β * t) * exp(-ξα * t)
    du[6] = n * exp(β * t) * exp(-ξα * t) * A
    du[7] = n * exp(β * t) * exp(ξα * t) * B
    du[8] = 0.0
end
function stoch_psa_solve_model(p::Tuple;tmax=10000.0,output=:all)
    α̂,ξα,ηα,β,rS,rR,dT,dD,s₀,r₀,p₀,σ,Tx = p
    # Solve the ODE
    sol = solve(ODEProblem(f!_stoch_psa_ode,[s₀,r₀,p₀,0.0,0.0,0.0,0.0,1],(0.0,tmax),p);verbose=false)
    # Process
    if output == :raw
        return sol
    elseif output == :psa
        return t -> sol(t)[3] / p₀
    else
        n = t -> sum(sol(t)[1:2])   # Cell population
        p̂ = t -> sol(t)[3] / p₀     # Normalised PSA
        A,B,C,D = [t -> sol(t)[i] for i = 4:7]
        Σ = (t,s) -> t ≤ s ? 
            exp(-β * (s + t)) * ηα^2 / (2 * ξα) / p₀^2 * (C(t) + A(t) * B(s) - D(t)) : 
            Σ(s,t)
        σ = t -> sqrt(max(0.0,Σ(t,t)))
        return n, p̂, σ, Σ
    end
end
stoch_psa_solve_model(q,Tx=default_Tx();kwargs...) = stoch_psa_solve_model(stoch_psa_get_pars(q,Tx);kwargs...)


#################################
## FITTING (EXACT LIKELIHOOD)
#################################

"""
    Get posteriors for a patient. By default, loads results from a saved JLD2 file.
"""
function stoch_psa_patient_fit(day, psa, Tx, iters=100000; guess = stoch_psa_prior, tryagain=5, kwargs...)

    # Setup loglikelihood
    function loglike(q)
        _, p̂, _, Σ = stoch_psa_solve_model(q,Tx;tmax=maximum(day),output=:all)
        μ = p̂.(day)
        σ = q[end]
        S = [Σ(d₁,d₂) for d₁ in day, d₂ in day] + σ^2 * I
        return try
            loglikelihood(MvNormal(μ,S), psa)
        catch e
            -Inf
        end
    end

    # Setup posterior
    logpost(q) = insupport(stoch_psa_prior,q) ?  loglike(q) + logpdf(stoch_psa_prior,q) : -Inf

    # Sample
    res = mcmc(guess, logpost, iters;names=stoch_psa_par_names, tryagain, kwargs...);

    # Check that we actually converged
    if any(rhat(res).nt.rhat .> 1.1)
        @warn("Chains did not converge!")
    end

    return res
    
end