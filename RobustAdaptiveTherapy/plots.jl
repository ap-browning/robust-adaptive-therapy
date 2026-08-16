"""
    plot_patient_fit(idxs;n=1000,kwargs...)

kwargs are passed to the `patient_fit` function.
"""
function plot_patient_fit(idx,res;n=1000,model=:strobl,show_prediction=false,kwargs...)
    
    # Load data, inference results
    day, psa, Tx = patient_data(idx)
    q = get_map(res)

    # p₀ varies between stochastic and deterministic models
    p₀_idx = findfirst(res.name_map.parameters .== :p₀)

    ## (a) : PSA

        # Function to get solution for each sample (PSA)
        func = q -> 
            t -> q[p₀_idx] * solve_model(q,Tx;tmax=maximum(day),output=:psa,model)(t)

        # Get quantiles, MAP
        T = range(0.0,maximum(day),201)
        l,u = get_quantiles(func,T,res,n)
        m = model != :stochastic ? func(q).(T) :    # Deterministic
            [func(q).(T) for q = sample(res,5)]     # Stochastic

        if show_prediction

            # Function to get a prediction for each sample
            func_pred = q -> 
                t -> rand(psa_noise_model(t, q)) + q[p₀_idx] * 
                    solve_model(q,Tx;tmax=maximum(day),output=:psa,model).(t)

            # Get quantiles for prediction interval
            lp,up = get_quantiles(func_pred, T, res, n; func_broadcast=false)

            # Plot...
            plt1 = plot(T,up,frange=lp,c=:orange,α=0.2,lw=0.0,label="95% PI")
            plot!(T,u,frange=l,c=:black,α=0.2,lw=0.0,label="95% CrI")
        else
            plt1 = plot(T,u,frange=l,c=:black,α=0.2,lw=0.0,label="95% CrI")
        end

        if model != :stochastic
            plot!(plt1,T,m,lw=2.0,c=:black,label="MAP")
        else
            plot!(plt1,T,m,lw=2.0,α=0.5,c=:black,label="")
        end
        scatter!(plt1,day, psa, c=:red, msw=0.0, label="Data")

    ## (b) : Treatment
    plt2 = plot(Tx, xlim=extrema(day),frange=0.0,lw=0.0,c=:blue,label="",yticks=[],ylabel="Tx")

    ## (c) : Resistant proportion

        # Function to get resistant proportion for each sample
        func = q -> solve_model(q,Tx;tmax=maximum(day),output=:all,model)[3]

        # Get quantiles, MAP
        l,u = get_quantiles(func,T,res,n)
        m = model != :stochastic ? func(q).(T) :    # Deterministic
            [func(q).(T) for q = sample(res,5)]     # Stochastic

        # Plot...
        plt3 = plot(T,u,frange=l,c=:black,α=0.2,lw=0.0,label="95% CrI")
        if model != :stochastic
            plot!(plt3,T,m,lw=2.0,c=:black,label="MAP")
        else
            plot!(plt3,T,m,lw=2.0,α=0.5,c=:black,label="")
        end

    ## (d): (Stochastic only) Cost of resistance
    if model == :stochastic

        # Function to obtain cost of resistance
        func = q -> 
            t -> solve_model(q,Tx;tmax=maximum(day),output=:raw,model)(t)[3]

        # Get quantiles
        l,u = get_quantiles(func,T,res,n)
        m = [func(q).(T) for q = sample(res,5)]

        # Plot...
        plt4 = plot(T,u,frange=l,c=:black,α=0.2,lw=0.0,label="95% CrI")
        plot!(plt4,T,m,lw=2.0,α=0.5,c=:black,label="")
        plot!(plt4,ylim=:auto, ylabel="Resist. Cost")

    end

    ## Figure...
    plot!(plt1,ylim=(0.0,1.2),ylabel="Norm. PSA")
    plot!(plt2,ylim=(0.0,1.0))
    plot!(plt3,ylim=(0.0,1.0),ylabel="Prop. Resist.")

    if model != :stochastic
        return plot(plt1,plt2,plt3,layout=@layout([a; b{0.1h}; c]),link=:x, box=:on,xlabel="Time [d]",plot_title="Patient $idx")
    else
        return plot(plt1,plt2,plt3,plt4,layout=@layout([a; b{0.1h}; c; d]),
            link=:x, box=:on,xlabel="Time [d]",plot_title="Patient $idx",size=(600,600))
    end
end
plot_patient_fit(idx;n=1000,kwargs...) = plot_patient_fit(idx,patient_fit(idx;kwargs...);n,kwargs...)


#################################
## USEFUL PLOTTING FUNCTIONS
#################################
"""
    get_quantiles(func,T,pars)
"""
function get_quantiles(func,T,pars;q=[0.025,0.975],func_broadcast=true)
    X = zeros(length(pars),length(T))
    for (i,par) = enumerate(pars)
        f = func(par[:])
        if func_broadcast
            X[i,:] = f.(T)
        else
            X[i,:] = f(T)
        end
    end
    ([quantile(x, qᵢ) for x in eachcol(X)] for qᵢ in q)
end
get_quantiles(func,T,C::Chains,n=1000;kwargs...) = get_quantiles(func,T,sample(C,n);kwargs...)


#################################
## OTHER PLOTS
#################################

"""
    density2d(x,y,...)

Create 2D kernel density plot.
"""
@userplot density2d
@recipe function f(kc::density2d; levels=10, clip=((-3.0, 3.0), (-3.0, 3.0)), z_clip = nothing)
    x,y = kc.args

    x = vec(x)
    y = vec(y)

    k = KernelDensity.kde((x, y))
    z = k.density / maximum(k.density)
    if !isnothing(z_clip)
        z[z .< z_clip * maximum(z)] .= NaN
    end

    legend --> false

    @series begin
        seriestype := contourf
        colorbar := false
        (collect(k.x), collect(k.y), z')
    end

end


function confellipse!(plt,x,y;level=0.95,remove_outliers=true,kwargs...)
    idx = trues(size(x))
    if remove_outliers
        for z = (x,y)
            q1,q2 = quantile(z)
            lb,ub = q1 - 1.5 * (q2 - q1), q2 + 1.5 * (q2 - q1)
            idx[@. !(lb < z < ub)] .= 0.0
        end
    end
    μ = [mean(x[idx]), mean(y[idx])]
    Σ = cov([x[idx] y[idx]])
    covellipse!(plt,μ,Σ;n_std = sqrt(quantile(Chisq(2),level)),kwargs...)
end
confellipse(x,y;kwargs...) = confellipse!(plot(), x, y; kwargs...)