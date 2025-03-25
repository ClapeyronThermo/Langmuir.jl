"""
    `ObrienMyers(M, K₀, σ, E)`

    ObrienMyers <: IsothermModel

`ObrienMyers(M, K₀, σ, E)` represents the O'Brien & Myers isotherm model.

## Inputs
⋅
- `M::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀::T`: Affinity parameter at T → ∞, `[Pa⁻¹]`
- `σ::T`: Variance in energy distribution, `[unitless]` #Not too sure...
- `E::T`: Adsorption energy, `[J⋅mol⁻¹]`

More on σ²: #Not too sure...
σ² is essentially the same σ² as is given in the expression of a normal distirbution. Basically, a larger σ² will correspond to a larger range of values taken by teh adsorption energies

## Description

The O'Brien & Myers equation is given by:

n = M⋅(K⋅p⋅(1+K⋅p)⁻¹+σ²⋅K⋅p⋅(1-K⋅p)⋅(2⋅(1+K⋅p)³)⁻¹)

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀*exp(-E / (R * T))

Where:
- `R` is the universal gas constant, `[J⋅mol⁻¹⋅K⁻¹]`,
- `T` is the temperature, `[K]`.
"""


@with_metadata struct ObrienMyers{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "affinity parameter") #Using Inf cause trouble in bboxoptimize
    (σ::T, (-Inf, Inf), "variance in energy distribution")
    (E::T, (-Inf, 0.0), "energy parameter")
end

function sp_res(model::ObrienMyers, p, T) #done TVC
    M = model.M
    K₀ = model.K₀
    σ = model.σ
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))
    _2 = convert(eltype(model), 2)
    return M * (log1p(K * p) + σ^(_2) * K * p / (_2 * (_1 + K * p)^(_2)))
end

function loading(model::ObrienMyers, p, T) #done TVC
    M = model.M
    K₀ = model.K₀
    σ = model.σ
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))
    _2 = convert(eltype(model), 2)
    _3 = convert(eltype(model), 3)
    return M * (K *p / (_1 + K*p) + σ^(_2)* K * p * (_1 - K * p) / (_2 * (_1 + K * p)^(_3))) #again not sure about _x instead of x...
end

#optimizations for ObrienMyers, not necessary, but improve performance | Done TVC
henry_coefficient(model::ObrienMyers, T) = model.M * model.K₀ * exp(-E / (Rgas(model) * T)) * (one(eltype(model)) + 0.5 * modoel.σ^(convert(eltype(model), 2)))
saturated_loading(model::AntiLangmuir, T) = model.M