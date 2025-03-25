"""
    `AntiLangmuir(A₀, B, E)`

    AntiLangmuir <: IsothermModel

`AntiLangmuir(A₀, B, E)` represents the Anti_Langmuir isotherm model.

## Inputs
⋅
- `A₀::T`: Henry coefficient, `[mol⋅kg⁻¹⋅Pa⁻¹]`
- `B::T`: Affinity parameter at T → ∞, `[Pa⁻¹]`
- `E::T`: Adsorption energy, `[J⋅mol⁻¹]`

## Description

The Anti-Langmuir equation is given by:

n = A⋅p/(1-B⋅p)

The adsorption energy E is related to the Henry coefficient A₀ by the equation:

A = A₀*exp(-E / (R * T))

Where:
- `R` is the universal gas constant, `[J⋅mol⁻¹⋅K⁻¹]`,
- `T` is the temperature, `[K]`.

"""


@with_metadata struct AntiLangmuir{T} <: IsothermModel{T}
    (A₀ ::T, (0.0, Inf), "Henry coefficient")
    (B::T, (0.0, Inf), "affinity parameter") #Using Inf cause trouble in bboxoptimize
    (E::T, (-Inf, 0.0), "energy parameter")
end

function sp_res(model::AntiLangmuir, p, T)
    A₀= model.A₀
    B = model.B
    E = model.E
    A = A₀*exp(-E/(Rgas(model)*T))
    return -A/B*log1p(-B*p) 
end

function loading(model::AntiLangmuir, p, T)
    A₀= model.A₀
    B = model.B
    E = model.E
    A = A₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))
    return A *p / (_1 - B*p)
end

#optimizations for Anti-Langmuir, not necessary, but improve performance
henry_coefficient(model::AntiLangmuir, T) = model.A₀*exp(-E/(Rgas(model)*T))
saturated_loading(model::AntiLangmuir, T) = prevfloat(Inf) #as temperature increases, so does the loading in all cases, thus the saturated loading rends to infinity
#saturated_loading(model::AntiLangmuir, T) = prevfloat(loading(model,Inf,Inf)) (not sure...)

export Anti_Langmuir