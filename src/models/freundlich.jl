@with_metadata struct Freundlich{T} <: IsothermModel{T}
    (K₀::T, (0.0, Inf), "Affinity parameter")
    (f₀::T, (0.0, Inf), "Surface heterogeneity parameter at T → ∞")
    (β::T, (0.0, Inf), "Surface heterogeneity coefficient")
    (E::T, (-Inf, 0.0), "Energy parameter")
end

function sp_res(model::Freundlich, p, T)
    K₀, f₀, β, E = model.K₀, model.f₀, model.β, model.E
    f = f₀ + β/T
    K = K₀*exp(-E/(Rgas(model)*T))
    return K*p^f/f
end

function pressure_impl(model::Freundlich, Π, T,::typeof(sp_res), approx)
    K₀, f₀, β, E = model.K₀, model.f₀, model.β, model.E
    f = f₀ + β/T
    K = K₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))
    v = _1/f
    return (Π/(K*v))^v
end

function x0_guess_fit(::Type{T}, data::AdsIsoTData) where T <: Freundlich
    l, p = data.l, data.p
    #l = M*p^f
    #log(l) = log(M) + f*log(p)
    logp, logl = log.(p), log.(l)
    _1 = one.(p)
    logK, f = hcat(_1,logp) \ logl
    K = exp(logK)
    return T(K, f, one(T), zero(T))
end

export Freundlich
