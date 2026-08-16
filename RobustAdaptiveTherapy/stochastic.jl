#################################
## PARAMETER HANDLING
#################################

"""
   `get_pars`` for the stochastic model
"""
function stoch_get_pars(q,Tx=default_Tx())
    KT,Kξ,Kη,rS, cR, d̂T, dD, n₀, fR, p₀, σ, τ = q
    rR = rS * cR
    dT = d̂T * rS
    s₀ = (1 - fR) * n₀
    r₀ = fR * n₀
    return (KT,Kξ,Kη,rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx)
end

"""
    `prior` for the stochastic model
"""
prior_KT = Truncated(Normal(1.0,0.5),0.0,Inf)
prior_Kξ = log(2) * InverseGamma(5,1/100)
prior_Kη = LogUniform(1e-4,1e-2)
stoch_prior = Product([[
    prior_KT, prior_Kξ, prior_Kη]; prior.v])

stoch_par_names = ["KT";"Kξ";"Kη"; par_names]

"""
    Convert from deterministic to stochastic parameters
"""
function stoch_convert_qode(qode)
    [mean(stoch_prior)[1:3]; qode]
end

#################################
## SDE MODEL
#################################

function f!_stoch(du,u,p,t)
    s,r,K,drug = u
    KT,Kξ,Kη,rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = p
    if isa(Tx,Function)
        drug = Tx(t)
    end
    du[1] = rS * s * (1 - (s + r) / K) * (1 - dD * drug) - dT * s
    du[2] = rR * r * (1 - (s + r) / K) - dT * r
    du[3] = -Kξ * (K - KT)
    du[4] = 0.0
end
function g!_stoch(du,u,p,t)
    du .= 0.0
    KT,Kξ,Kη,rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = p
    du[3] = Kη
end


"""
    stoch_solve_model(p;tmax=10000.0,output=:psa,model=:strobl)

Solve the model. The tuple p contains both the parameters and information about the treatment schedule. For example:
    `p = get_pars(q)`
which uses (by default) the default treatment strategy.


    solve_model(q,Tx=default_Tx();kwargs...)

Alternative calling sequence.
"""
function stoch_solve_model(p::Tuple;tmax=10000.0,output=:psa,dt=1.0)
    KT,Kξ,Kη,rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = p
    # Solve the ODE
    if isa(Tx,Function)
        # Tx is given as a function of time
        sol = solve(SDEProblem(f!_stoch,g!_stoch,[s₀,r₀,1.0,1.0],(0.0,tmax),p),EM(),dt=dt;verbose=false)
    else
        # Tx is given as a callback (i.e., strategy)
        sol = solve(SDEProblem(f!_stoch,g!_stoch,[s₀,r₀,1.0,1.0],(0.0,tmax),p),EM(),dt=dt;verbose=false,callback=Tx)
    end
    psa = t -> sum(sol(t)[1:2]) / (s₀ + r₀)
    drug = t -> sol(t)[4] .> 0.5
    rprp = t -> sol(t)[2] / sum(sol(t)[1:2])
    if output == :raw
        return sol
    elseif output == :psa
        return psa
    else
        return psa, drug, rprp
    end
end
stoch_solve_model(q,Tx=default_Tx();kwargs...) = stoch_solve_model(stoch_get_pars(q,Tx);kwargs...)

#################################
## INFERENCE
#################################

"""
    Create particle filter approximation to the likelihood
"""
function stoch_create_pm_loglike(day, psa, Tx)

    function pm_loglike(q;nparticle=200)

        # Parameters
        p = stoch_get_pars(q,Tx)
        KT,Kξ,Kη,rS,rR,dT,dD,s₀,r₀,p₀,σ,τ,Tx = p
        u₀ = [s₀,r₀,1,1]

        # Setup problem
        prob = SDEProblem(f!_stoch, g!_stoch, u₀, (0.0,maximum(day)), p)
         
        # Step from t₀ to t₁ with u(t₀) = u₀
        function step(u0;tspan=(0.0,1.0))
            sol = solve(remake(prob;u0,tspan),EM(),dt=1.0,saveat=tspan)
            return sol.u[end]
        end

        # Observation model (independent / t = 0)
        function y(x)
            Normal(p₀ * sum(x[1:2]) / (s₀ + r₀), σ)
        end

        # Observation model (dependent)
        ρ = δt -> τ == 0.0 ? 0.0 : exp(-abs(δt) * log(2) / τ)
        p̂ = x -> p₀ * sum(x[1:2]) / (s₀ + r₀)
        function y(xnew,xold,tnew,told,pold)
            p̂old = p̂(xold)
            p̂new = p̂(xnew)
            ρδt = ρ(tnew - told)
            Normal(p̂new + ρδt * (pold - p̂old), σ * sqrt(1 - ρδt^2))
        end

        # Initialise random particles
        x = [u₀ for _ = 1:nparticle]
        xold = deepcopy(x)

        # Initialise loglike
        loglike = log(mean(pdf.(y.(x),psa[1])))

        # Step through times...
        for i = eachindex(day)[2:end]

            # Save old particles
            @. xold = x

            # Step all particles forward
            x = step.(x;tspan=day[i-1:i])

            # Reweight particles
            Wt = pdf.(y.(x,xold,day[i],day[i-1],psa[i-1]),psa[i])
            if all(Wt .== 0.0) || any(isinf.(Wt)) || any(isnan.(Wt))
                # Particle degeneracy
                return -Inf
            end
            W = Weights(Wt)
            
            # Resample
            x = sample(x,W,nparticle)

            # Update loglikelihood
            loglike += log(mean(W))

        end

        return loglike

    end

    return pm_loglike

end

"""
    Stochastic version of `patient_fit`
"""
function stoch_patient_fit(idx,iters=100000;nparticle=100,load_saved=true,save_result=true,tryagain=0,kwargs...)

    # Check that the patient exists
    idx ∈ patient_list || error("Patient doesn't exist!")

    if load_saved

        fname = "$(@__DIR__)/fits/stochastic/fit$idx.jld2"
        if isfile(fname)
            @load fname res
            return res
        else
            @warn("File for this patient doesn't exist, running inference again...")
        end

    end

    # Load data
    day, psa, Tx = patient_data(idx)

    # Fix for patient 71, who has a duplicate measurement on one day (identical)
    if idx == 71
        idup = findfirst(day[1:end-1] .== day[2:end])
        idup_rem = setdiff(eachindex(day),idup)
        day = day[idup_rem]
        psa = psa[idup_rem]
    end

    # For some patients, ignore initial measurement (likely incorrectly recorded)
    if idx ∈ [1,7,11,15,18,95,101]
        fit_idx = eachindex(day)[2:end]
        day = day[fit_idx]
        psa = psa[fit_idx]
    end

    # Setup loglikelihood (particle filter)
    loglike = stoch_create_pm_loglike(day, psa, Tx)

    # Setup posterior
    logpost(q) = insupport(stoch_prior,q) ?  loglike(q;nparticle) + logpdf(stoch_prior,q) : -Inf

    # Initial guess (try to load deterministic results and use MAP)
    qsde = stoch_convert_qode(get_map(patient_fit(idx;model=:strobl)))

    # Sample
    res = mcmc(qsde, logpost, iters;names=stoch_par_names,tryagain);

    # Check that we actually converged
    if any(rhat(res).nt.rhat .> 1.1)
        @warn("Chains did not converge for patient $idx. Results will still be saved for stochastic model!")
    end

    # Save
    if save_result
        @save "$(@__DIR__)/fits/stochastic/fit$idx.jld2" res
        fig = plot_patient_fit(idx;model=:stochastic)
        savefig(fig,"$(@__DIR__)/fits/stochastic/fit$idx.png")
    end
    return res

end