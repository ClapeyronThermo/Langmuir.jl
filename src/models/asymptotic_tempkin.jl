"""
    `AsymptoticTempkin(M, K₀, θ, E)`

    AsymptoticTempkin <: IsothermModel

`AsymptoticTempkin(M, K₀, θ, E)` represents the asymptotic Tempkin isotherm model.

## Inputs

- `M::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `θ::T`: Correction factor, `[unitless]`
- `E::T`: Adsorption energy, `[J⋅mol⁻¹]`

## Description

The Asymptotic Tempkin equation is given by:

n = M ⋅ (K ⋅ p ÷ (1 + K ⋅ p) + θ ⋅ (K ⋅ p ÷ (1 + K ⋅ p))² ⋅ (K ⋅ p ÷ (1 + K ⋅ p) - 1)

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀ ⋅ exp(-E ÷ (R ⋅ T))

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct AsymptoticTempkin{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading") 
    (K₀::T, (0.0, Inf), "affinity parameter") #Using Inf cause trouble in bboxoptimize
    (θ::T, (0.0, 1.0), "correction factor") 
    (E::T, (-Inf, 0.0), "energy parameter") 
end

function sp_res(model::AsymptoticTempkin, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    θ = model.θ
    _1 = one(eltype(model))
    _2 = two(eltype(model))
    return M * (log1p(K * p) - _1 / _2  * θ * (K * p / (_1 + K * p))^(_2)) #not sure if 0.5 would work instead of (_1/_2)...
end


function loading(model::AsymptoticTempkin, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    θ = model.θ
    _1 = one(eltype(model))
    _2 = two(eltype(model))

    return M * (K * p / (_1 + K * p) + θ * (K * p / (_1 + K * p))^(_2) * (K * p / (_1 + K * p) - _1))
end

#optimizations for AsymptoticTempkin, not necessary, but improve performance
henry_coefficient(model::AsymptoticTempkin, T) = model.K₀*exp(-E/(Rgas(model)*T)) * model.M
saturated_loading(model::AsymptoticTempkin, T) = model.M 