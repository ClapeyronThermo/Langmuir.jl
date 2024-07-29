function fit(data::AdsIsoTData, ::Type{M}, x0 = to_vec(x0_guess_fit(M, data))) where M <: IsothermModel
    p,l,T = getfield(data,1),getfield(data,2),getfield(data,3)
    function ℓ(θ)
		# construct model
		model = from_vec(M, θ)
		loss = zero(eltype(model))
        p,l = data.p, data.l
		for (pᵢ, nᵢ, Tᵢ) in zip(p,l,T)
			# predicted loading
            n̂ᵢ = loading(model, pᵢ, Tᵢ)
			# increment loss.
			loss += (nᵢ - n̂ᵢ)^2
		end
		return loss
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
    result = optimize(ℓ,x0,NLSolvers.LineSearch(NLSolvers.LBFGS()))
    return from_vec(M, x_sol(result)),x_minimum(result)
end