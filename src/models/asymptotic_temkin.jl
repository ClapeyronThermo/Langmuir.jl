"""
    `AsymptoticTemkin(M, K₀, θ, E)`

    AsymptoticTemkin <: IsothermModel

`AsymptoticTemkin(M, K₀, θ, E)` represents the asymptotic Temkin isotherm model.

## Inputs

- `M::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `θ::T`: Strength of adsorbate-adsorbate interactions, `[unitless]`
- `E::T`: Adsorption energy, `[J⋅mol⁻¹]`

## Description

The value of θ can be both negative and positive. 
In the former case, this would correspond to attraction between adsorbed particle, and in the latter to repulsion.

The Asymptotic Temkin equation is given by:

n = M ⋅ (K ⋅ p ÷ (1 + K ⋅ p) + θ ⋅ (K ⋅ p ÷ (1 + K ⋅ p))² ⋅ (K ⋅ p ÷ (1 + K ⋅ p) - 1)

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀ ⋅ exp(-E ÷ (R ⋅ T))

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct AsymptoticTemkin{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading") 
    (K₀::T, (0.0, Inf), "affinity parameter") #Using Inf cause trouble in bboxoptimize
    (θ::T, (-Inf, Inf), "strength of adsorbate-adsorbate interactions") 
    (E::T, (-Inf, 0.0), "energy parameter") 
end

function sp_res(model::AsymptoticTemkin, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    θ = model.θ
    _1 = one(eltype(model))
    _2 = convert(eltype(model), 2)
    return M * (log1p(K * p) - _1 / _2  * θ * (K * p / (_1 + K * p))^(_2)) 
end


function loading(model::AsymptoticTemkin, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    θ = model.θ
    _1 = one(eltype(model))
    _2 = convert(eltype(model), 2)

    return M * (K * p / (_1 + K * p) + θ * (K * p / (_1 + K * p))^(_2) * (K * p / (_1 + K * p) - _1))
end

#optimizations for AsymptoticTemkin, not necessary, but improve performance
henry_coefficient(model::AsymptoticTemkin, T) = model.K₀*exp(-E/(Rgas(model)*T)) * model.M
saturated_loading(model::AsymptoticTemkin, T) = model.M 

export AsymptoticTemkin