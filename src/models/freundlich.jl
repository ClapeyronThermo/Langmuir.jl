"""
    `Freundlich(K₀, f₀, β, E)`

     Freundlich <: IsothermModel


## Inputs

- `K₀::T`: Affinity parameter at T → ∞, `[1/Pa]`
- `f₀::T`: Surface heterogeneity parameter at T → ∞, `[-]`
- `β::T`: Surface heterogeneity coefficient `[K]`
- `E::T`: Adsorption energy, `[J/mol]`

## Description

The Freundlich isotherm is given by:

n = K_f × pᶠ

The affinity parameter `K_f` is a temperature dependent and can be linked to adsorption energy `E` by:

K_f = K₀ × exp(-E / (RT))

The exponent f is also temperature dependent and can be expressed as: 

f = f₀ - β/T

Where:
- `R` is the gas constant,
- `T` is the temperature.
"""
@with_metadata struct Freundlich{T} <: IsothermModel{T}
    (K₀::T, (0.0, Inf), "Affinity parameter")
    (f₀::T, (0.0, Inf), "Surface heterogeneity parameter at T → ∞")
    (β::T, (0.0, Inf), "Surface heterogeneity coefficient")
    (E::T, (-Inf, 0.0), "Energy parameter")
end

function sp_res(model::Freundlich, p, T)
    K₀, f₀, β, E = model.K₀, model.f₀, model.β, model.E
    f = f₀ - β/T
    K = K₀*exp(-E/(Rgas(model)*T))
    return K*p^f/f
end

function pressure_impl(model::Freundlich, Π, T,::typeof(sp_res))
    K₀, f₀, β, E = model.K₀, model.f₀, model.β, model.E
    f = f₀ - β/T
    K = K₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))
    v = _1/f
    return (Π/(K*v))^v
end

function x0_guess_fit(::Type{T}, data::AdsIsoTData) where T <: Freundlich

    Ts, l_p = split_data_by_temperature(data)

    logKs = Vector{eltype(Ts)}(undef, length(l_p))
    fs = Vector{eltype(Ts)}(undef, length(l_p))
    _1_ = ones(eltype(Ts), length(first(first(l_p))))


    _1 = one(eltype(Ts))
    _0 = zero(eltype(Ts))
    _1s = ones(eltype(Ts), length(Ts))

    for i in eachindex(l_p)
        l_i, p_i = l_p[i]
        logl_i, logp_i = log.(l_i), log.(p_i)
        logKs[i], fs[i] = hcat(_1_, logp_i) \ logl_i
    end
    
    #l = K*p^f
    #log(l) = log(K) + f*log(p)


    if length(l_p) > 1
        logK0, E = hcat(_1s, _1./ (Rgas(T).*Ts)) \ logKs
        f0, β = hcat(_1s, -_1./Ts) \ fs
        K0 = exp(logK0)
    else
        K0 = first(exp(logKs))
        f0 = _1
        β = _0
        E = _1
    end


    return T(K0, f0, β, -E)
end

export Freundlich
