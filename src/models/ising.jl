"""
    `IsingS1(Mᵢ,Kᵢ,Eᵢ,Kₒ,Eₒ)`

    IsingS1 <: IsothermModel

`IsingS1(Mᵢ,Kᵢ,Eᵢ,Kₒ,Eₒ)` represents the single site Ising isotherm model.

## Inputs

- `Mᵢ::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `Kᵢ::T`: Affinity parameter XX at T → ∞, `[Pa⁻¹]`
- `Eᵢ::T`: Adsorption energy XX, `[J⋅mol⁻¹]`
- `Kₒ::T`: Affinity parameter XX at T → ∞, `[Pa⁻¹]`
- `Eₒ::T`: Adsorption energy XX, `[J⋅mol⁻¹]`


## Description

The IsingS1 equation is given by:

n = 

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀*exp(-E / (R * T))

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct IsingS1{T} <: IsothermModel{T}
    (Mᵢ ::T, (0.0, Inf), "saturation loading")
    (Kᵢ::T, (0.0, Inf), "affinity parameter I") #Using Inf cause trouble in bboxoptimize
    (Eᵢ::T, (-Inf, 0.0), "energy parameter I")
    (Kₒ::T, (0.0, Inf), "affinity parameter O") #Using Inf cause trouble in bboxoptimize
    (Eₒ::T, (-Inf, 0.0), "energy parameter O")
end

'XXX'

function loading(model::LangmuirS1, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))

    return M * K *p / (_1 + K*p)
end

#optimizations for LangmuirS1, not necessary, but improve performance
henry_coefficient(model::LangmuirS1, T) = model.M*model.K₀*exp(-model.E/(Rgas(model)*T))
saturated_loading(model::LangmuirS1, T) = model.M #Some depend on T, some don't
pressure_impl(model::LangmuirS1, Π, T,::typeof(sp_res), approx) = expm1(Π/model.M)/(model.K₀*exp(-model.E/(Rgas(model)*T)))

