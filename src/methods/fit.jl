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
function IsothermFittingProblem(IsothermModel::Type{M}, loading_data::AdsIsoTData{TT}, loss::L = abs2; fittable::Union{Nothing,AbstractVector{Bool}}=nothing) where {M <: IsothermModel, TT, L}
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
    seed::Int = 314159265358
end

Base.@kwdef struct NewtonIsothermFittingSolver <: IsothermFittingSolver
    logspace::Bool = false
end

Base.@kwdef struct NLSolversIsothermFittingSolver{T} <: IsothermFittingSolver
    method::T = NLSolvers.BFGS()
    verbose::Bool = false
    logspace::Bool = false
    ftol::Float64 = 1e-8
    xtol::Float64 = 1e-8
end

#= function heterogeneity_penalty(M::Union{Freundlich, LangmuirFreundlich, Sips, RedlichPeterson}, T)
       
    return ifelse(M.f₀ - M.β/T < 0.0, 50.0, zero(eltype(M))) 

end =#

#from (0,1) to (lb,ub)
function auto_interp_eval(_x,_lb,_ub,logspace = false)
    x,lb,ub = promote(_x,_lb,_ub)
    
    x = clamp(x,zero(x),one(x))

    if isinf(lb) && isinf(ub)
        y = logspace ? atanh(2*x - 1) : tanpi(x - 0.5)
        return y
    end
    #translation
    w = abs(ub) > abs(lb) ? max(zero(lb),-2*lb) : min(zero(ub),-2*ub)
    ubt,lbt = ub + w, lb + w

    if isinf(lb) && !isinf(ub)
        #=
        from -Inf to ub
        at x = 0 -> y = -Inf
        at x = 1 -> y = ub
        =#
        is_log = (logspace || log2(abs(ub)) < -16) && !iszero(ub)
        x̄ = is_log ? 1/(1 - log2(x)) : x
        y = -1/x̄ - ubt*x̄ + 1 + 2*ub + w

    elseif !isinf(lb) && isinf(ub)
        #=
        from lb to Inf
        at x = 0 -> y = lb
        at x = 1 -> y = Inf
        =#
        is_log = (logspace || log2(abs(lb)) < -16) && !iszero(lb)
        x̄ = is_log ? 1 - 1/(1 - log2(x)) : 1 - x
        y = 1/x̄ - lbt*x̄ -1 + 2*lb + w

    else
        ratio = ubt/lbt
        maybe_log_negative = ubt < 0 && ratio < 0.01
        maybe_log_positive = lbt > 0 && ratio > 100
        if logspace || maybe_log_positive || maybe_log_negative
            if iszero(ub) #[ub,0]
                lny = (1 - x)*log2(-lb + 1)
                y = 1 - exp2(lny)
            elseif iszero(lb) #[0,ub]
                lny = x*log2(ub + 1)
                y = exp2(lny) - 1
            else #[lb,ub]
                y = lbt*exp2(x*log2(ratio)) - w
            end
        else
            y = x*ub + (1 - x)*lb
        end
    end
    return y
end


#from (lb,ub) to (0,1)
function auto_interp_inv(_y, _lb, _ub, logspace = false)
    y,lb,ub = promote(_y,_lb,_ub)
    y = clamp(y,lb,ub)
    @assert lb < ub "invalid bounds, lb = $lb is equal or greater than ub = $ub"

    if isinf(lb) && isinf(ub)
        y = logspace ? 0.5*(tanh(y) + 1) : atan(y)/π + 0.5
        return y
    end

    #translation
    w = abs(ub) > abs(lb) ? max(zero(lb),-2*lb) : min(zero(ub),-2*ub)
    ubt,lbt = ub + w, lb + w
    yt = y + w
    if isinf(lb) && !isinf(ub) 
        #y = -1/x̄  - ubt*x̄   +  1 + 2*ub + w
        #0 = -1/x̄  - ubt*x̄   +  1 + 2*ub + w - y
        #0 =   -1  - ubt*x̄^2 + (1 + 2*ub + w - y)*x̄
        a = -ubt
        b = (1 + 2*ub + w - y)
        c = -1
        if !iszero(a)
            disc2 = b^2 - 4*a*c
            disc = sqrt(disc2)
            x̄₁ = (-b + disc) / (2*a)
            x̄₂ = (-b - disc) / (2*a)
            x̄ = (0 < x̄₁ < 1) ? x̄₁ : x̄₂
        else
            x̄ = -c/b
        end
        is_log = (logspace || log2(abs(ub)) < -16) && !iszero(ub)
        x = is_log ? exp2(1 - 1/x̄) : x̄

    elseif !isinf(lb) && isinf(ub)
        #y = 1/x̄ - lbt*x̄      - 1 + 2*lb + w
        #0 = 1/x̄ - lbt*x̄      - 1 + 2*lb + w - y
        #0 = 1   - lbt*x̄^2 + (- 1 + 2*lb + w - y)*x̄
        a = -lbt
        b = (- 1 + 2*lb + w - y)
        c = 1
        if !iszero(a)
            disc2 = b^2 - 4*a*c
            disc = sqrt(disc2)
            x̄₁ = (-b + disc) / (2*a)
            x̄₂ = (-b - disc) / (2*a)
            x̄ = (0 < x̄₁ < 1) ? x̄₁ : x̄₂
        else
            x̄ = -c/b
        end
        u = 1 - x̄
        is_log = (logspace || log2(abs(lb)) < -16) && !iszero(lb)
        x = is_log ? exp2(1 - 1/u) : u

    else
        ratio = ubt/lbt
        maybe_log_negative = ubt < 0 && ratio < 0.01
        maybe_log_positive = lbt > 0 && ratio > 100
        if logspace || maybe_log_positive || maybe_log_negative
            if iszero(ub) #[ub,0]
                u = log2(1 - y)/log2(-lb + 1)
                x = 1 - u
            elseif iszero(lb) #[0,ub]
                x = log2(y + 1)/log2(ub + 1)
            else
                x = log2(yt/lbt)/log2(ratio)
            end
        else
            # y = x*ub + (1-x)*lb  →  x = (y - lb)/(ub - lb)
            x = (y - lb) / (ub - lb)
        end
    end
    return x
end

#translate from [lb,ub] to [-Inf,Inf]
function unconstrain(_y,_lb,_ub,logspace = true)
    y,lb,ub = promote(_y,_lb,_ub)
    isinf(lb) && isinf(ub) && (return y)
    x = auto_interp_inv(y,lb,ub,logspace)
    return auto_interp_eval(x,-Inf,Inf,logspace)
end

#translate from [-Inf,Inf] to [lb,ub]
function constrain(_ȳ,_lb,_ub,logspace = true)
    ȳ,lb,ub = promote(_ȳ,_lb,_ub)
    isinf(lb) && isinf(ub) && (return ȳ)
    x = auto_interp_inv(ȳ,-Inf,Inf,logspace)
    return auto_interp_eval(x,lb,ub,logspace)
end

function CommonSolve.solve(prob::IsothermFittingProblem{M, L, DL, DC, X, LB, UB, F},
alg::DEIsothermFittingSolver) where {M, L, DL, DC, X, LB, UB, F}
    
    Ðₗ = prob.LoadingData
    Ðₕ = prob.CalorimetricData

    p = pressure(Ðₗ)
    l = loading(Ðₗ)
    T = temperature(Ðₗ)
    σ² = variance(Ðₗ)
    _0 = zero(Base.promote_eltype(p,l,T,σ²))
    _1 = one(_0)
    n_max,i_max = findmax(l)
    n_max += 3*sqrt(σ²[i_max]) #maximum bound
    ℓ(θ) = let prob = prob, p = p, l = l, T = T, is_logspace = alg.logspace, n_max = n_max,ub = prob.ub, lb = prob.lb

        _θ = copy(θ)
        _θ .= auto_interp_eval.(θ,lb,ub,is_logspace)

        # Reconstruct full model from fittable parameters
        model = from_vec_fittable(prob.IsothermModel, _θ, prob.model_template, prob.fittable)

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

    x0 = auto_interp_inv.(prob.x0,prob.lb,prob.ub,alg.logspace)
    seed0 = BlackBoxOptim.Random.seed!()
    BlackBoxOptim.Random.seed!(alg.seed)
    
    result = BlackBoxOptim.bboptimize(ℓ, x0; 
                                    SearchRange = [(_0, _1) for i in eachindex(x0)],
                                    PopulationSize = alg.population_size,
                                    MaxTime = alg.time_limit,
                                    MaxSteps = alg.max_steps,
                                    TraceMode = ifelse(alg.verbose, :verbose, :silent))

    BlackBoxOptim.Random.seed!(seed0)

    x_best = BlackBoxOptim.best_candidate(result)

    loss_opt_M = BlackBoxOptim.best_fitness(result)

    θ_best = similar(x_best)
    θ_best .= auto_interp_eval.(x_best,prob.lb,prob.ub,alg.logspace)

    # Reconstruct full model from fitted fittable parameters
    return loss_opt_M, from_vec_fittable(prob.IsothermModel, θ_best, prob.model_template, prob.fittable)
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
    method = alg.method
    
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
    alg) where {M, L, DL, DC, X, LB, UB, F}
    
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
             solver=DEIsothermFittingSolver()) where M <: IsothermModel
    
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