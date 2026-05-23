module LangmuirMetaheuristicsExt
using Langmuir
using Metaheuristics
using Langmuir.CommonSolve

function CommonSolve.solve(prob::IsothermFittingProblem{M, L, DL, DC, X, LB, UB, F},
alg::Metaheuristics.AbstractAlgorithm) where {M, L, DL, DC, X, LB, UB, F}
    
    Ðₗ = prob.LoadingData
    Ðₕ = prob.CalorimetricData

    p = Langmuir.pressure(Ðₗ)
    l = Langmuir.loading(Ðₗ)
    T = Langmuir.temperature(Ðₗ)
    σ² = Langmuir.variance(Ðₗ)
    _0 = zero(Base.promote_eltype(p,l,T,σ²))
    _1 = one(_0)
    n_max,i_max = findmax(l)
    n_max += 3*sqrt(σ²[i_max]) #maximum bound
    ℓ(θ) = let prob = prob, p = p, l = l, T = T, n_max = n_max,ub = prob.ub, lb = prob.lb

        _θ = copy(θ)
        _θ .= Langmuir.auto_interp_eval.(θ,lb,ub)

        # Reconstruct full model from fittable parameters
        model = Langmuir.from_vec_fittable(prob.IsothermModel, _θ, prob.model_template, prob.fittable)

        ℓr = zero(eltype(model))
        
        for (pᵢ, nᵢ, Tᵢ, σ²ᵢ) in zip(p, l, T, σ²)
            n̂ᵢ = loading(model, pᵢ, Tᵢ) #Predicted loading
            ℓrᵢ = iszero(σ²ᵢ) ? prob.loss(nᵢ - n̂ᵢ) : prob.loss(nᵢ - n̂ᵢ)/σ²ᵢ #if zero variance, just use the loss
            #!isfinite(ℓrᵢ) && (ℓrᵢ += 1e100) #if not finite, add a big number
            #n̂ᵢ < 0 && (ℓrᵢ *= exp(100*abs(n̂ᵢ))) #if loading is negative, add a penalty multiplier, proportional to the negativity
            #n̂ᵢ > n_max && (ℓrᵢ *= exp(abs(n_max - n̂ᵢ))) #if loading is greater than maximum loading, also add penalty
            ℓr += ℓrᵢ
        end

        return ℓr

    end

    x0 = Langmuir.auto_interp_inv.(prob.x0,prob.lb,prob.ub)
    Metaheuristics.set_user_solutions!(alg, x0, ℓ);
    bounds = Metaheuristics.boxconstraints(lb = zeros(length(x0)), ub =ones(length(x0)))
    result = Metaheuristics.optimize(ℓ,bounds,alg)

    x_best = Metaheuristics.minimizer(result)

    loss_opt_M = Metaheuristics.minimum(result)

    θ_best = similar(x_best)
    θ_best .= Langmuir.auto_interp_eval.(x_best,prob.lb,prob.ub)

    # Reconstruct full model from fitted fittable parameters
    return loss_opt_M, Langmuir.from_vec_fittable(prob.IsothermModel, θ_best, prob.model_template, prob.fittable)
end

end