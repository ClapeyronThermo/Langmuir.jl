#general multisite langmuir.

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


#optimizations for Langmuir, not necessary, but improve performance
henry_coefficient(model::Langmuir, T) = model.M*model.K₀*exp(-model.E/(Rgas(model)*T))
saturated_loading(model::Langmuir, T) = model.M #Some depend on T, some don't
sp_res_pressure_impl(model::Langmuir, Π, T) = expm1(Π/model.M)/(model.K₀*exp(-model.E/(Rgas(model)*T)))


"""
    DualSiteLangmuir <: IsothermModel

    DualSiteLangmuir(M1, K₀1, E1, M2, K₀2, E2)

 DualSiteLangmuir(M1, K₀1, E1, M2, K₀2, E2) represents the DualSite Langmuir isotherm model.

## Inputs
- `M1`::T: maximum loading capacity of the first adsorbent site, `[mol/kg]`
- `K₀1`::T: equilibrium constant at zero coverage for the first site, `[1/Pa]`
- `E1`::T: adsorption energy for the first site, `[J/mol]`
- `M2`::T: maximum loading capacity of the second adsorbent site, `[mol/kg]`
- `K₀2`::T: equilibrium constant at zero coverage for the second site, `[1/Pa]`
- `E2`::T: adsorption energy for the second site, `[J/mol]`

## Description

The DualSite Langmuir equation is given by:

n = M1 * log(1 + K1 * p) + M2 * log(1 + K2 * p)

where:
- n is the loading of the adsorbate on the adsorbent,
- M1 and M2 are the maximum loading capacities of the first and second adsorbent sites, respectively,
- K1 and K2 are the equilibrium constants at zero coverage for the first and second sites, respectively,
- p is the pressure of the gas.

The adsorption energies E1 and E2 are related to the equilibrium constants K₀1 and K₀2 by the equations:

K1 = exp(-E1 / (R * T))
K2 = exp(-E2 / (R * T))
    
where:
- R is the gas constant,
- T is the temperature,

"""
struct DualSiteLangmuir{T} <: IsothermModel{T}
    M1::T
    K₀1::T
    E1::T
    M2::T
    K₀2::T
    E2::T
end

function sp_res(model::DualSiteLangmuir, p, T)
    M1, K₀1, E1, M2, K₀2, E2 = model.M1, model.K₀1, model.E1, model.M2, model.K₀2, model.E2
    RT = Rgas(model)*T
    K1 = K₀1*exp(-E1/RT)
    K2 = K₀2*exp(-E2/RT)
    return M1*log1p(K1*p) + M2*log1p(K2*p)
end

#optimizations for DualSiteLangmuir, not necessary, but improve performance

function henry_coefficient(model::DualSiteLangmuir, T) 
    RT = Rgas(model)*T    
    return model.M1*model.K₀1*exp(-model.E1/RT) + model.M2*model.K₀2*exp(-model.E2/RT)
end

saturated_loading(model::DualSiteLangmuir, T) = model.M1 + model.M2

export Langmuir, DualSiteLangmuir, henry_coefficient