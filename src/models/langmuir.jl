"""
    `Langmuir(M, K₀, E)`

    Langmuir <: IsothermModel

Langmuir(M, K₀, E) represents the Langmuir isotherm model, which describes the adsorption of a gas on a solid surface.

## Inputs

- `M`::T: maximum loading capacity of the adsorbent, `[mol/kg]`
- `K₀`::T: equilibrium constant at zero coverage, `[1/Pa]`
- `E`::T: adsorption energy, `[J/mol]`

## Description

The Langmuir equation is given by:

n = (M * K₀ * p) / (1 + K₀ * p)

where:
- n is the loading of the adsorbate on the adsorbent,
- M is the maximum loading capacity of the adsorbent,
- K₀ is the equilibrium constant at zero coverage,
- p is the pressure of the gas.

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K₀ = exp(-E / (R * T))

where:
- R is the gas constant,
- T is the temperature.

"""
struct Langmuir{T} <: IsothermModel{T}
    M::T
    K₀::T
    E::T
end

function sp_res(model::Langmuir, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    return M*log1p(K*p)
end


function loading(model::Langmuir, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))

    return M * K *p / (_1 + K*p)
end

#optimizations for Langmuir, not necessary, but improve performance
henry_coefficient(model::Langmuir, T) = model.M*model.K₀*exp(-model.E/(Rgas(model)*T))
saturated_loading(model::Langmuir, T) = model.M #Some depend on T, some don't
pressure_impl(model::Langmuir, Π, T,::typeof(sp_res), approx) = expm1(Π/model.M)/(model.K₀*exp(-model.E/(Rgas(model)*T)))

#TODO: include effects of temperature. at the moment, the fit procedure ignores temperature dependence.
#probably requires separating the models by temperature and linearizing K to obtain T-dependence.
function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: Langmuir
    # use first two data points to get the slope


    #l = M*k*p/(1 + k*p)
    #l*(1 + k*p) = M*k*p
    #l + l*k*p = M*k*p
    #M*k*p - l*k*p = l
    #p*(Mk) - l*p(k) = l
    
    l,p = data.l,data.p
    MK,K = hcat(p,-l .* p)\l
    M = MK/K
    return Langmuir(M,K,zero(K))
end

export Langmuir
