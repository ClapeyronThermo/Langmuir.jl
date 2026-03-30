abstract type IsothermFittingSolver end

struct IsothermFittingProblem{M <: IsothermModel, TT, DL, DC, L, X, LB, UB, F}
    IsothermModel::Type{M}
    LoadingData::AdsIsoTData{TT}
    CalorimetricData::DC
    loss::L
    x0::X
    lb::LB
    ub::UB
    fittable::F
    model_template::M
end


# Constructor with all arguments
IsothermFittingProblem(IsothermModel::Type{M}, loading_data::AdsIsoTData{TT}, calorimetric_data::DC, loss::L, x0::X, lb::LB, ub::UB, fittable::F, model_template::M) where {M <: IsothermModel, TT, DC, L, X, LB, UB, F} = 
IsothermFittingProblem{M, TT, typeof(loading_data), DC, L, X, LB, UB, F}(IsothermModel, loading_data, calorimetric_data, loss, x0, lb, ub, fittable, model_template)

# Simplified constructor
function IsothermFittingProblem(IsothermModel::Type{M}, loading_data::AdsIsoTData{TT}, loss::L; fittable::Union{Nothing,AbstractVector{Bool}}=nothing) where {M <: IsothermModel, TT, L}
    model_template = x0_guess_fit(IsothermModel, loading_data)
    
    # Default: all parameters are fittable
    if fittable === nothing
        fittable = trues(model_length(IsothermModel))
    end
    
    # Validate fittable length
    length(fittable) == model_length(IsothermModel) || throw(ArgumentError("fittable vector length must match number of model parameters"))
    
    # Extract fittable parameters
    x0 = to_vec_fittable(model_template, fittable)
    full_lb = isotherm_lower_bound(eltype(loading_data), IsothermModel)
    full_ub = isotherm_upper_bound(eltype(loading_data), IsothermModel)
    
    # Get bounds for only fittable parameters
    fittable_indices = findall(fittable)
    lb = [full_lb[i] for i in fittable_indices]
    ub = [full_ub[i] for i in fittable_indices]
    
    return IsothermFittingProblem(IsothermModel, loading_data, nothing, loss, x0, lb, ub, fittable, model_template)
end


Base.@kwdef struct DEIsothermFittingSolver <: IsothermFittingSolver
    max_steps::Int = 2e4
    population_size::Int = 500
    time_limit::Float64 = 20.0
    verbose::Bool = false
    logspace::Bool = true
end

Base.@kwdef struct NewtonIsothermFittingSolver <: IsothermFittingSolver
    logspace::Bool = false
end

Base.@kwdef struct NLSolversIsothermFittingSolver <: IsothermFittingSolver
    max_iter::Int = 1000
    verbose::Bool = false
    logspace::Bool = false
    ftol::Float64 = 1e-8
    xtol::Float64 = 1e-8
end

#= function heterogeneity_penalty(M::Union{Freundlich, LangmuirFreundlich, Sips, RedlichPeterson}, T)
       
    return ifelse(M.f₀ - M.β/T < 0.0, 50.0, zero(eltype(M))) 

end =#

function CommonSolve.solve(prob::IsothermFittingProblem{M, L, DL, DC, X, LB, UB, F},
alg::DEIsothermFittingSolver) where {M, L, DL, DC, X, LB, UB, F}
    
    Ðₗ = prob.LoadingData
    Ðₕ = prob.CalorimetricData

    p = pressure(Ðₗ)
    l = loading(Ðₗ)
    T = temperature(Ðₗ)
    σ² = variance(Ðₗ)
    lb_sign = sign.(nextfloat.(prob.lb)) #Assumes that parameters are either > 0 or < 0.

     ℓ(θ) = let prob = prob, lb_sign = lb_sign, p = p, l = l, T = T, is_logspace = alg.logspace

        _θ = copy(θ)

        if is_logspace
            _θ .= lb_sign .* exp.(θ)
        end
        
        # Reconstruct full model from fittable parameters
        model = from_vec_fittable(prob.IsothermModel, _θ, prob.model_template, prob.fittable)

        ℓr = zero(eltype(model))

        for (pᵢ, nᵢ, Tᵢ, σ²ᵢ) in zip(p, l, T, σ²)

            n̂ᵢ = loading(model, pᵢ, Tᵢ) #Predicted loading

            if isnan(n̂ᵢ)
                n̂ᵢ = -one(nᵢ)*nᵢ
            end

			ℓr += prob.loss(nᵢ - n̂ᵢ)/σ²ᵢ
        end

        return ℓr

    end


    if alg.logspace
        x0 = log.(sign.(prob.x0) .* prob.x0)
        lb = log.(sign.(nextfloat.(prob.lb)) .* nextfloat.(prob.lb))
        ub = log.(sign.(prevfloat.(prob.ub)) .* prevfloat.(prob.ub))
    
        # Ensure lb and ub are mutable arrays
        lb = collect(lb)
        ub = collect(ub)
    
        # Swap lb and ub in logspace if the parameter is between -Inf and 0
        for i in eachindex(lb)
            if lb[i] > ub[i]
                lb[i], ub[i] = ub[i], lb[i]
            end
        end
        
    else

        lb = nextfloat.(prob.lb)
        ub = prevfloat.(prob.ub)
        x0 = prob.x0 
    end


    result = BlackBoxOptim.bboptimize(ℓ, x0; 
    SearchRange = [(lb[i], ub[i]) for i in eachindex(lb)],
    PopulationSize = alg.population_size,
    MaxTime = alg.time_limit,
    MaxSteps = alg.max_steps,
    TraceMode = ifelse(alg.verbose, :verbose, :silent))

    opt_M = BlackBoxOptim.best_candidate(result)

    loss_opt_M = BlackBoxOptim.best_fitness(result)

    opt_θ = similar(opt_M)

    if alg.logspace
        opt_θ .= lb_sign .* exp.(opt_M)
    else
        opt_θ = opt_M
    end

    # Reconstruct full model from fitted fittable parameters
    return loss_opt_M, from_vec_fittable(prob.IsothermModel, opt_θ, prob.model_template, prob.fittable)

end

function CommonSolve.solve(solver::NewtonIsothermFittingSolver)
        nothing #TODO
end

function CommonSolve.solve(prob::IsothermFittingProblem{M, L, DL, DC, X, LB, UB, F},
alg::NLSolversIsothermFittingSolver) where {M, L, DL, DC, X, LB, UB, F}
    
    Ðₗ = prob.LoadingData
    p = pressure(Ðₗ)
    l = loading(Ðₗ)
    T = temperature(Ðₗ)
    σ² = variance(Ðₗ)
    
    # Box constraints for NLSolvers
    lb = collect(prob.lb)
    ub = collect(prob.ub)
    
    # Objective function (sum of squared errors)
    function loss_fn(θ)
        # Reconstruct full model from fittable parameters
        model = from_vec_fittable(prob.IsothermModel, θ, prob.model_template, prob.fittable)
        
        ℓr = zero(eltype(model))
        for (pᵢ, nᵢ, Tᵢ, σ²ᵢ) in zip(p, l, T, σ²)
            n̂ᵢ = loading(model, pᵢ, Tᵢ)
            if isnan(n̂ᵢ)
                n̂ᵢ = -one(nᵢ) * nᵢ
            end
            ℓr += prob.loss(nᵢ - n̂ᵢ) / σ²ᵢ
        end
        return ℓr
    end
    
    # Use NLSolvers with box constraints
    x0 = copy(prob.x0)
    
    # Set up optimization options
    options = NLSolvers.OptimizationOptions(
        maxiter = alg.max_iter,
        f_abstol = alg.ftol,
        x_abstol = alg.xtol,
        show_trace = alg.verbose
    )
    
    # Use L-BFGS method - for bounded optimization, use through TrustRegion/ActiveBox
    method = NLSolvers.BFGS()
    
    # Create objective with automatic differentiation
    obj = ADScalarObjective(loss_fn, x0)
    
    # Create optimization problem with bounds (inplace=true required)
    prob_opt = NLSolvers.OptimizationProblem(obj, (lb, ub); inplace = true)
    
    # Solve the optimization problem starting from x0
    result = NLSolvers.solve(prob_opt, x0, method, options)
    
    opt_θ = result.info.solution
    loss_opt = result.info.fx
    
    # Reconstruct full model from fitted fittable parameters
    return result, loss_opt, from_vec_fittable(prob.IsothermModel, opt_θ, prob.model_template, prob.fittable)
end

function fit(prob::IsothermFittingProblem{M, L, DL, DC, X, LB, UB, F},
    alg::IsothermFittingSolver) where {M, L, DL, DC, X, LB, UB, F}
    
    return solve(prob, alg)

end

"""
    fit(::Type{M}, data::AdsIsoTData; fittable=nothing, loss=abs2, solver=DEIsothermFittingSolver()) where M <: IsothermModel

Fit an isotherm model to adsorption data.

## Arguments
- `M`: The isotherm model type to fit
- `data`: Adsorption isotherm data
- `fittable`: Optional `Vector{Bool}` indicating which parameters to fit (default: all parameters are fittable)
- `loss`: Loss function (default: `abs2`)
- `solver`: Fitting solver (default: `DEIsothermFittingSolver()`)
  - `DEIsothermFittingSolver()`: Global stochastic optimizer (Differential Evolution)
  - `NLSolversIsothermFittingSolver()`: Local deterministic optimizer (L-BFGS)

## Returns
- `(loss_value, fitted_model)`: Tuple of the final loss value and the fitted model

## Example
```julia
# Fit all parameters with default stochastic solver
loss, model = fit(Freundlich, data)

# Fit only K₀, f₀, and E (keep β fixed at its initial value)
loss, model = fit(Freundlich, data, fittable=[true, true, false, true])

# Use deterministic L-BFGS solver
loss, model = fit(Freundlich, data, solver=NLSolversIsothermFittingSolver())
```
"""
function fit(::Type{M}, data::AdsIsoTData; 
             fittable::Union{Nothing,AbstractVector{Bool}}=nothing,
             loss=abs2,
             solver::IsothermFittingSolver=DEIsothermFittingSolver()) where M <: IsothermModel
    
    prob = IsothermFittingProblem(M, data, loss; fittable=fittable)
    return fit(prob, solver)
end

#= function fit(data::AdsIsoTData{TT}, ::Type{M}, loss = abs2, x0 = to_vec(x0_guess_fit(M, data))) where {M <: IsothermModel,TT}
    p = pressure(data)
    l = loading(data)
    T = temperature(data)
    function ℓ(θ)
		# construct model
		model = from_vec(M, θ)
		ℓr = zero(eltype(model))
		for (pᵢ, nᵢ, Tᵢ) in zip(p, l, T)
			# predicted loading
            n̂ᵢ = loading(model, pᵢ, Tᵢ)
			# increment loss.
            if isnan(n̂ᵢ)
                n̂ᵢ = -one(nᵢ)*nᵢ
            end
			ℓr += loss(nᵢ - n̂ᵢ)
        end
		return ℓr
	end
    #=
    some notes:
    on fit, we depend on a function defined over θ,p,T ∈ ℝ
    while loading(model(θ),p,T) could be defined over θ,p,T ∈ ℝ,
    this is not true for sp_res(model(θ),p,T).

    for example. for LangmuirS1, loading is defined over θ,p,T ∈ ℝ
    but sp_res has the restriction K*p > -1
    if loading is calculated via AD, we have to keep an eye on those situations.
    =#
    lb = similar(x0)
    ub = similar(x0)
    lb .= isotherm_lower_bound(eltype(x0),M)
    ub .= isotherm_upper_bound(eltype(x0),M)
    result = optimize(ℓ,x0,NLSolvers.ActiveBox(factorize=NLSolvers.positive_factorize),bounds = (lb,ub))
    return from_vec(M, x_sol(result)),x_minimum(result)
end =#


export IsothermFittingProblem, DEIsothermFittingSolver, NLSolversIsothermFittingSolver, fit