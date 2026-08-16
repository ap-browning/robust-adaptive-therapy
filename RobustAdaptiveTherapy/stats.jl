#################################
## MCMC
#################################

"""
    Run MCMC to sample from a distribution proportional to `fdensity`.
"""
function mcmc(q,fdensity,iters;
    names=["p$i" for i in 1:length(q)],
    algorithm=:aswam,nchain=4,tryagain=0,L=2,map_restart=false,kwargs...)

    # Adaptive MCMC
    res = Array{Any}(undef,nchain)
    @threads for i = 1:nchain
        if isa(q, Distribution)
            res[i] = adaptive_rwm(rand(q), fdensity, iters; algorithm, L, kwargs...)
        else
            res[i] = adaptive_rwm(q, fdensity, iters; algorithm, L, kwargs...)
        end
    end

    # Create chain
    C = Chains(permutedims(cat([r.X for r in res]...,dims=3),[2,1,3]),names,
        evidence=cat([r.D[1] for r in res]...,dims=2)')

    # Check that we actually moved!
    if length(unique(C.logevidence)) == 1
        # Try again...
        display("Trying again!")
        return mcmc(q,fdensity,iters;names,algorithm,nchain,tryagain,kwargs...)
    else
        if any(rhat(C).nt.rhat .> 1.1)
            if tryagain > 0
                display("Trying again! ($tryagain)")
                restart_location = map_restart ? get_map(C) : q
                return mcmc(restart_location,fdensity,iters;names,algorithm,nchain,tryagain=tryagain-1,kwargs...)
            else
                @warn("Chains may not have converged, max(R̂) > 1.1.")
            end
        end

        return C
    end

end

"""
    rhatval(C::Chains)
"""
rhatval(C::Chains) = rhat(C).nt.rhat
isconverged(C::Chains,thresh=1.1) = all(rhatval(C) .< thresh)

"""
    get_map(C::Chains)

Obtain the posterior mode.
"""
function get_map(C::Chains)
    idx = findmax(C.logevidence)[2]
    C.value[idx[2], :, idx[1]]
end

"""
    Convert chains to a matrix
"""
function Base.Matrix(C::Chains)
    vcat([C.value.data[:,:,i] for i = 1:size(C.value.data)[3]]...)
end

"""
    Resample from a set of chains
"""
function StatsBase.sample(C::Chains, n::Int)
    M = Matrix(C)
    rand([eachrow(M)...],n)
end


#################################
## INFERENCE - PRIOR
#################################

prior_rS = Truncated(Normal(0.027, 0.027 * 0.5),0.0,Inf)
prior_r̂R = Uniform(0.0,1.0)
prior_d̂T = Uniform(0.0,1.0)
prior_dD = Uniform(0.0,2.0)
prior_n₀ = Uniform(0.1,1.0)
prior_fR = Uniform(0.001,0.5)
prior_p₀ = Uniform(0.0,2.0)
prior_σ = Uniform(0.001,1.0)
prior_τ = Exponential(10.0)

prior = Product([prior_rS, prior_r̂R, prior_d̂T, prior_dD, prior_n₀, prior_fR, prior_p₀, prior_σ, prior_τ])
 
par_names = ["rS", "r̂R", "d̂T", "dD", "n₀", "fR", "p₀", "σ", "τ"]

#################################
## FITTING
#################################



"""
    Get posteriors for a patient. By default, loads results from a saved JLD2 file.
"""
function patient_fit(day, psa, Tx, iters=100000; guess = prior, model=:strobl, tryagain=5, kwargs...)

    if model == :stochastic
        return stoch_patient_fit(idx,iters;load_saved,save_result,tryagain,kwargs...)
    elseif model == :stochastic_psa
        return stoch_psa_patient_fit(day, psa, Tx, iters; tryagain, kwargs...)
    end

    # Setup loglikelihood
    function loglike(q)
        # Initial PSA level
        p₀ = q[end-2]
        # Solve model
        sol = solve_model(q,Tx;tmax=maximum(day),output=:psa,model)
        # Likelihood
        logpdf(psa_noise_model(day,q) + p₀ * sol.(day),psa)
    end

    # Setup posterior
    logpost(q,day,psa,Tx) = insupport(prior,q) ?  loglike(q) + logpdf(prior,q) : -Inf

    # Sample
    res = mcmc(guess, q -> logpost(q, day, psa, Tx),iters;names=par_names,tryagain, kwargs...);

    # Check that we actually converged
    if any(rhat(res).nt.rhat .> 1.1)
        @warn("Chains did not converge!")
    end

    return res
    
end

function patient_fit(idx,iters=100000;load_saved=true,save_result=true, tryagain=5, model=:strobl, skip_save = true, kwargs...)

    if model == :stochastic
        return stoch_patient_fit(idx,iters;load_saved,save_result,tryagain,kwargs...)
    elseif model == :stochastic2
        return stoch2_patient_fit(idx,iters;load_saved,save_result,tryagain,kwargs...)
    end

    # Check that the patient exists
    idx ∈ patient_list || error("Patient doesn't exist!")

    if load_saved

        fname = "$(@__DIR__)/fits/$(String(model))/fit$idx.jld2"
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
    if idx ∈ [1,7,11,14,15,18,28,95,101]
        fit_idx = eachindex(day)[2:end]
        day = day[fit_idx]
        psa = psa[fit_idx]
    end

    res = patient_fit(day, psa, Tx, iters; model, tryagain, kwargs...)
    
    # Check that we actually converged
    if any(rhat(res).nt.rhat .> 1.1) && (skip_save == true)
        @warn("Chains did not converge for patient $idx. Results will not be saved!")
        return res
    end

    if any(rhat(res).nt.rhat .> 1.1)
        saveloc = "$(@__DIR__)/fits/$(String(model))/unobtained"
    else
        saveloc = "$(@__DIR__)/fits/$(String(model))"
    end

    # Save
    if save_result
        @save "$saveloc/fit$idx.jld2" res
        fig = plot_patient_fit(idx,res)
        savefig(fig,"$saveloc/fit$idx.png")
    end
    return res

end


"""
    patient_fit_list()

Get all patients that have posterior distributions saved, that have 
"""
function patient_fit_list(;mindata=true,model=:strobl,progressonly=false,all_patients=false)
    dir = "$(@__DIR__)/fits/$(String(model))"
    files = filter(x -> endswith(x,"jld2"),readdir(dir))
    lst = sort([parse(Int64,split(x,".")[1][4:end]) for x in files])
    if progressonly
        if model == :stochastic
            lst = intersect(lst, patient_fit_list(model=:strobl,progressonly=true) )
        else
            lst = lst[model_patient_map_progress.(lst)]
        end
    end
    if mindata
        return intersect(lst, patient_list_mindata)
    else
        return lst
    end
end

"""
    get_map(idx::Number)

Get the best-fit parameter set for a given patient
"""
get_map(idx::Number;model=:strobl) = get_map(patient_fit(idx;model))

"""
    sample_joint_posterior(n;mindata=true,storedonly=true)
"""
function sample_joint_posterior(n=1;mindata=true,storedonly=true,model=:strobl,kwargs...)
    lst = storedonly ? patient_fit_list(;mindata=false,model,kwargs...) : patient_list
    lst = mindata ? intersect(lst,patient_list_mindata) : lst
    # 1: Sample patients (minimise loading of data for efficiency)
    samp_pts_all = rand(lst,n)
    samp_pts_idx = unique(samp_pts_all)
    samp_pts_num = [count(samp_pts_all .== idx) for idx = samp_pts_idx]

    # 2: Look through patients and sample the required number
    vcat([sample(patient_fit(samp_pts_idx[i];model),samp_pts_num[i]) for i = eachindex(samp_pts_idx)]...)
end
