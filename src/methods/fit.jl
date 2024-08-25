abstract type IsothermFittingSolver end

struct IsothermFittingProblem{M <: IsothermModel, TT, DL, DC, L, X, LB, UB}
    IsothermModel::Type{M}
    LoadingData::AdsIsoTData{TT}
    CalorimetricData::DC
    loss::L
    x0::X
    lb::LB
    ub::UB
end


# Constructor with all arguments
IsothermFittingProblem(IsothermModel::Type{M}, loading_data::AdsIsoTData{TT}, calorimetric_data::DC, loss::L, x0::X, lb::LB, ub::UB) where {M <: IsothermModel, TT, DC, L, X, LB, UB} = 
IsothermFittingProblem{M, TT, typeof(loading_data), DC, L, X, LB, UB}(IsothermModel, loading_data, calorimetric_data, loss, x0, lb, ub)

# Simplified constructor
IsothermFittingProblem(IsothermModel::Type{M}, loading_data::AdsIsoTData{TT}, loss::L) where {M <: IsothermModel, TT, L} = 
    IsothermFittingProblem(IsothermModel, loading_data, nothing, loss,
        to_vec(x0_guess_fit(IsothermModel, loading_data)),
        isotherm_lower_bound(eltype(loading_data), IsothermModel),
        isotherm_upper_bound(eltype(loading_data), IsothermModel))


Base.@kwdef struct DEIsothermFittingSolver <: IsothermFittingSolver
    max_steps::Int = 2e4
    population_size::Int = 50
    time_limit::Float64 = Inf
    verbose::Bool = false
    logspace::Bool = true
end

Base.@kwdef struct NewtonIsothermFittingSolver <: IsothermFittingSolver
    logspace::Bool = false
end

function CommonSolve.solve(prob::IsothermFittingProblem{M, L, DL, DC, X, LB, UB},
alg::DEIsothermFittingSolver) where {M, L, DL, DC, X, LB, UB}
    
    Ðₗ = prob.LoadingData
    Ðₕ = prob.CalorimetricData

    p = pressure(Ðₗ)
    l = loading(Ðₗ)
    T = temperature(Ðₗ)
    lb_sign = sign.(prob.lb) #Assumes that parameters are either > 0 or < 0.

    function ℓ(θ)

        _θ = copy(θ)

        if alg.logspace
            _θ .= lb_sign .* exp.(θ)
        end
        
        model = from_vec(prob.IsothermModel, _θ)

        ℓr = zero(eltype(model))

        for (pᵢ, nᵢ, Tᵢ) in zip(p, l, T)

            n̂ᵢ = loading(model, pᵢ, Tᵢ) #Predicted loading

            if isnan(n̂ᵢ)
                n̂ᵢ = -one(nᵢ)*nᵢ
            end

			ℓr += prob.loss(nᵢ - n̂ᵢ)
        end

        return ℓr/length(p)

    end


    if alg.logspace

        x0 = log.(sign.(prob.x0) .* prob.x0)
        lb = log.(sign.(prob.lb) .* prob.lb)
        ub = log.(sign.(prob.ub) .* prob.ub)
    
        # Ensure lb and ub are mutable arrays
        lb = collect(lb)
        ub = collect(ub)
    
        # Swap lb and ub in logspace if the parameter is between -Inf and 0
        for i in eachindex(lb)
            if prob.lb[i] < 0 && prob.ub[i] < 0
                lb[i], ub[i] = ub[i], lb[i]
            end
    end
        
        else
            lb = prob.lb
            ub = prob.ub
            x0 = prob.x0
        
    end


    result = BlackBoxOptim.bboptimize(ℓ, x0; 
    SearchRange = [(lb[i], ub[i]) for i in eachindex(lb)],
    PopulationSize = alg.population_size,
    MaxTime = alg.time_limit,
    TraceMode = ifelse(alg.verbose, :verbose, :silent))

    opt_M = BlackBoxOptim.best_candidate(result)

    loss_opt_M = BlackBoxOptim.best_fitness(result)

    opt_θ = similar(opt_M)

    if alg.logspace
        opt_θ .= lb_sign .* exp.(opt_M)
    else
        opt_θ = opt_M
    end

    return loss_opt_M, from_vec(prob.IsothermModel, opt_θ)

end

function CommonSolve.solve(solver::NewtonIsothermFittingSolver)
        nothing
end

function fit(prob::IsothermFittingProblem{M, L, DL, DC, X, LB, UB},
    alg::IsothermFittingSolver) where {M, L, DL, DC, X, LB, UB}
    
    return solve(prob, alg)

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

    for example. for langmuir, loading is defined over θ,p,T ∈ ℝ
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


export IsothermFittingProblem, DEIsothermFittingSolver, fit