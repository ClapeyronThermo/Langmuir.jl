
"""
    Toth <: IsothermModel

    Toth(M, K₀, E, f₀, β)

## Inputs

 - `M`::T: maximum loading capacity of the adsorbent, `[mol/kg]`
 - `K₀::T: equilibrium constant at zero coverage, `[1/Pa]`
 - `E`::T: adsorption energy, `[J/mol]`
 - `f₀`::T: Empirical parameter, `-`
 - `β`::T: Empirical parameter, `K`

## Description
K = K₀*exp(-E/(RT))
f = f₀ + β/T
nᵢ = M*K*P/(1 + (K*P)ᶠ)¹/ᶠ

"""
struct Toth{T} <: IsothermModel{T}
    M::T
    K₀::T
    E::T
    f₀::T
    β::T
end


function loading(model::Toth, p, T)
    M = model.M
    K = model.K₀*exp(-model.E/(Rgas(model)*T))
    f = model.f₀ + model.β/T
    return M*K*p/(1 + (K*p)^f)^(1/f)
end


henry_coefficient(model::Toth, T) = model.M*model.K₀*exp(-model.E/(Rgas(model)*T))
saturated_loading(model::Toth, T) = model.M #Some depend on T, some don't


export Toth, henry_coefficient