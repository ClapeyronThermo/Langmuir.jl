"""
    `RedlichPeterson(M, K₀, E, f)`

    RedlichPeterson <: IsothermModel

RedlichPeterson(M, K₀, E, f) represents the Redlich-Peterson isotherm model, which describes the adsorption of a gas on a solid surface.

## Inputs

- `M`::T: maximum loading capacity of the adsorbent, `[mol/kg]`
- `K₀`::T: equilibrium constant at zero coverage, `[1/Pa]`
- `E`::T: adsorption energy, `[J/mol]`
- `f₀`::T: parameter characterising the heterogeneity of the system  (no units)
- `β`::T: coefficient characterising the heterogeneity of the system  (K)

## Description

The RedlichPeterson equation is given by:

n = M * p / (1 + (K₀ * p)^f)

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀ × exp(-E / (R * T))

The exponent f is also temperature dependent and can be expressed as: 

f = f₀ - β/T

where:
- R is the gas constant,
- T is the temperature.
"""
@with_metadata struct RedlichPeterson{T} <: IsothermModel{T}
    (M::T , (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "affinity parameter") 
    (E::T, (-Inf, 0.0), "energy parameter")
    (f₀::T, (0.0, Inf), "surface heterogeneity parameter at T → ∞")
    (β::T, (0.0, Inf), "surface heterogeneity coefficient")
end

function loading(model::RedlichPeterson, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    f = model.f₀ - model.β/T
    Kpf = K*p^f
    _1 = one(eltype(model))
    return M*p/(_1 + Kpf)
end

#TODO: sp_res is defined, but depends on Hypergeometric functions.

henry_coefficient(model::RedlichPeterson, T) = model.M
#pressure_impl(model::RedlichPeterson, Π, T,::typeof(sp_res), approx) = expm1(Π/model.M)/(model.K₀*exp(-model.E/(Rgas(model)*T)))

function x0_guess_fit(::Type{T}, data::AdsIsoTData) where T <: RedlichPeterson
    langmuir_model = x0_guess_fit(LangmuirS1,data)
    M, K₀, E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E
    _0 = zero(M)
    _1 = one(M)
    return T(M, K₀, E, _1, _0)
end


export RedlichPeterson
