module LangmuirMetaheuristicsExt
using Langmuir
using Metaheuristics
using Langmuir.CommonSolve

function CommonSolve.solve(prob::IsothermFittingProblem{M, L, DL, DC, X, LB, UB, F},
alg::Metaheuristics.AbstractAlgorithm) where {M, L, DL, DC, X, LB, UB, F}
    data = prob.LoadingData
    _0 = zero(Base.promote_eltype(pressure(data),loading(data),temperature(data),variance(data),prob.x0,prob.lb,prob.ub))
    _1 = one(_0)
    #n_max,i_max = findmax(l)
    #n_max += 3*sqrt(σ²[i_max]) #maximum bound
    ℓ(θ) = let prob = prob
        Langmuir.isotherm_fitting_loss(prob,θ)
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