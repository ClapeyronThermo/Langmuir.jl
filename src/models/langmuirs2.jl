"""
    `LangmuirS2(M₁, K₀₁, E₁, M₂, K₀₂, E₂)`

    LangmuirS2 <: IsothermModel

`LangmuirS2(M₁, K₀₁, E₁, M₂, K₀₂, E₂)` represents the 2 site Langmuir isotherm model.

## Inputs

- `M₁::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₁::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₁::T`: Adsorption energy, `[J⋅mol⁻¹]`
- `M₂::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₂::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₂::T`: Adsorption energy, `[J⋅mol⁻¹]`

## Description

The LangmuirS2 equation is given by:

n = (M₁ ⋅ K₁ ⋅ p) ÷ (1 + K₁ ⋅ p) + (M₂ ⋅ K₂ ⋅ p) ÷ (1 + K₂ ⋅ p)

The adsorption energy Eᵢ is related to the equilibrium constant K₀ᵢ by the equation:      

Kᵢ = K₀ᵢ ⋅ exp(-Eᵢ ÷ (R ⋅ T))

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct LangmuirS2{T} <: IsothermModel{T}
    (M₁::T, (0.0, Inf), "saturation loading 1")
    (K₀₁::T, (0.0, Inf), "affinity parameter 1") #Using Inf cause trouble in bboxoptimize
    (E₁::T, (-Inf, 0.0), "energy parameter 1")
    (M₂::T, (0.0, Inf), "saturation loading 2")
    (K₀₂::T, (0.0, Inf), "affinity parameter 2") #Using Inf cause trouble in bboxoptimize
    (E₂::T, (-Inf, 0.0), "energy parameter 2")
end

function sp_res(model::LangmuirS2, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    return M₁*log1p(K₁*p)+M₂*log1p(K₂*p)
end


function loading(model::LangmuirS2, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    _1 = one(eltype(model))
    return M₁ * K₁ *p / (_1 + K₁*p)+ M₂ * K₂ *p / (_1 + K₂*p)
end


henry_coefficient(model::LangmuirS2, T) = model.M₁*model.K₀₁*exp(-model.E₁/(Rgas(model)*T))+model.M₂*model.K₀₂*exp(-model.E₂/(Rgas(model)*T))
saturated_loading(model::LangmuirS2, T) = model.M₁ + model.M₂

export LangmuirS2