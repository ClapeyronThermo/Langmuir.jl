"""
    `RedlichPeterson(M, K₀, E, f)`

    RedlichPeterson <: IsothermModel

RedlichPeterson(M, K₀, E, f) represents the Redlich-Peterson isotherm model, which describes the adsorption of a gas on a solid surface.

## Inputs

- `M`::T: maximum loading capacity of the adsorbent, `[mol/kg]`
- `K₀`::T: equilibrium constant at zero coverage, `[1/Pa]`
- `E`::T: adsorption energy, `[J/mol]`
- `f`::T: parameter characterising the heterogeneity of the system  (no units)

## Description

The RedlichPeterson equation is given by:

n = M * p / (1 + (K₀ * p)^f)

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K₀ = exp(-E / (R * T))

where:
- R is the gas constant,
- T is the temperature.

"""
struct RedlichPeterson{T} <: IsothermModel{T}
    M::T
    K₀::T
    E::T
    f::T
end

function loading(model::RedlichPeterson, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    Kpf = K*p^f
    return M*p/(1 + Kpf)
end

#TODO: sp_res is defined, but depends on Hypergeometric functions.

henry_coefficient(model::RedlichPeterson, T) = model.M
#pressure_impl(model::RedlichPeterson, Π, T,::typeof(sp_res), approx) = expm1(Π/model.M)/(model.K₀*exp(-model.E/(Rgas(model)*T)))

export RedlichPeterson
