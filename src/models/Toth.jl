
"""
    `Toth(M, K₀, E, f₀, β)`

    Toth <: IsothermModel


## Inputs

 - `M::T`: Saturation loading, `[mol/kg]`
 - `K₀::T`: Affinity parameter at T → ∞, `[Pa⁻¹]`
 - `E::T`: Adsorption energy, `[J/mol]`
 - `f₀::T`: Surface heterogeneity parameter at T → ∞, `[-]`
 - `β::T`: Surface heterogeneity coefficient, `[K]`

## Description

Toth isotherm model: 

n = M*K*P/(1 + (K*P)ᶠ)¹/ᶠ

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀ × exp(-E / (R * T))

The exponent f is also temperature dependent and can be expressed as: 

f = f₀ - β/T

where:
- R is the gas constant,
- T is the temperature.


"""
@with_metadata struct Toth{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "affinity parameter")
    (E::T, (-Inf, 0.0), "energy parameter")
    (f₀::T, (0.0, Inf), "surface heterogeneity parameter at T → ∞")
    (β::T, (-Inf, Inf), "surface heterogeneity coefficient")
end

function loading(model::Toth, p, T)
    M = model.M
    K = model.K₀*exp(-model.E/(Rgas(model)*T))
    f = model.f₀ - model.β/T
    Kpf = abs(K*p)^f #f has to be between 0 and 1
    _1 = one(eltype(p))
    return M*K*p/(_1 + Kpf)^(_1/f)
end

henry_coefficient(model::Toth, T) = model.M*model.K₀*exp(-model.E/(Rgas(model)*T))
saturated_loading(model::Toth, T) = model.M #Some depend on T, some don't

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: Toth
    langmuir_model = x0_guess_fit(LangmuirS1,data)
    M, K₀, E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E    
    _0 = nextfloat(zero(M))
    _1 = one(M)
    return T(M, K₀, E, _1, _0)
end

export Toth
