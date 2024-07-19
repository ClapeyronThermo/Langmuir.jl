#general multisite langmuir.

"""

Langmuir <: IsothermModel

Langmuir(M, K₀, E)

## Inputs

 - `M`::T: maximum loading capacity of the adsorbent, `[mol/kg]`
 - `K₀::T: equilibrium constant at zero coverage, `[1/Pa]`
 - `E`::T: adsorption energy, `[J/mol]`

## Description

A Langmuir isotherm model:

```math
q = M \frac{Kp}{1 + Kp}
K = K₀ \exp\left(\frac{-E}{RT}\right)
``` 

"""


struct Langmuir{T} <: IsothermModel{T}
    M::T
    K₀::T
    E::T
end

"""

sp_res(model::IsothermModel, p, T)

default units: `[mol/kg]`

Calculates the single component spreading pressure of the `model` given the bulk pressure `p` and temperature `T`:

## Returns

Returns a scalar of the same type as `p`.

Mathematically:

```math
```

"""


function sp_res(model::Langmuir, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    return M*log1p(K*p)
end


"""
henry_coefficient(model)

default units: `[mol/kg]`


Returns the single component spreading pressure of the `model` given the bulk pressure `p`, defined as:

```math
H = 

```


"""
#optimizations for Langmuir, not necessary, but improve performance
henry_coefficient(model::Langmuir, T) = model.M*model.K₀*exp(-model.E/(Rgas(model)*T))
saturated_loading(model::Langmuir, T) = model.M #Some depend on T, some don't
sp_res_pressure_impl(model::Langmuir, Π, T) = expm1(Π/model.M)/(model.K₀*exp(-model.E/(Rgas(model)*T)))

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

export Langmuir, DualSiteLangmuir