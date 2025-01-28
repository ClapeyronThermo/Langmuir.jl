"""
    `LangmuirS3(M₁, K₀₁, E₁, M₂, K₀₂, E₂, M₃, K₀₃, E₃)`

    LangmuirS3 <: IsothermModel

`LangmuirS3(M₁, K₀₁, E₁, M₂, K₀₂, E₂, M₃, K₀₃, E₃)` represents the single 3 Langmuir isotherm model.

## Inputs

- `M₁::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₁::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₁::T`: Adsorption energy, `[J⋅mol⁻¹]`
- `M₂::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₂::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₂::T`: Adsorption energy, `[J⋅mol⁻¹]`
- `M₃::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₃::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₃::T`: Adsorption energy, `[J⋅mol⁻¹]`

## Description

The LangmuirS3 equation is given by:

n = (M₁ ⋅ K₁ ⋅ p) ÷ (1 + K₁ ⋅ p) + (M₂ ⋅ K₂ ⋅ p) ÷ (1 + K₂ ⋅ p) + (M₃ ⋅ K₃ ⋅ p) ÷ (1 + K₃ ⋅ p)

The adsorption energy Eᵢ is related to the equilibrium constant K₀ᵢ by the equation:      

Kᵢ = K₀ᵢ ⋅ exp(-Eᵢ ÷ (R ⋅ T))

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct LangmuirS3{T} <: IsothermModel{T}
    (M₁::T, (0.0, Inf), "saturation loading 1")
    (K₀₁::T, (0.0, Inf), "affinity parameter 1") 
    (E₁::T, (-Inf, 0.0), "energy parameter 1")
    (M₂::T, (0.0, Inf), "saturation loading 2")
    (K₀₂::T, (0.0, Inf), "affinity parameter 2") 
    (E₂::T, (-Inf, 0.0), "energy parameter 2")
    (M₃::T, (0.0, Inf), "saturation loading 3")
    (K₀₃::T, (0.0, Inf), "affinity parameter 3") 
    (E₃::T, (-Inf, 0.0), "energy parameter 3")
end

function sp_res(model::LangmuirS3, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    M₃ = model.M₃
    K₀₃ = model.K₀₃
    E₃ = model.E₃
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    K₃ = K₀₃*exp(-E₃/(Rgas(model)*T))
    return M₁*log1p(K₁*p)+M₂*log1p(K₂*p)+M₃*log1p(K₃*p)
end


function loading(model::LangmuirS3, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    M₃ = model.M₃
    K₀₃ = model.K₀₃
    E₃ = model.E₃
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    K₃ = K₀₃*exp(-E₃/(Rgas(model)*T))
    _1 = one(eltype(model))
    return M₁ * K₁ *p / (_1 + K₁*p)+ M₂ * K₂ *p / (_1 + K₂*p)+ M₃ * K₃ *p / (_1 + K₃*p)
end

henry_coefficient(model::LangmuirS3, T) = model.M₁*model.K₀₁*exp(-model.E₁/(Rgas(model)*T))+model.M₂*model.K₀₂*exp(-model.E₂/(Rgas(model)*T))+model.M₃*model.K₀₃*exp(-model.E₃/(Rgas(model)*T))
saturated_loading(model::LangmuirS3, T) = model.M₁ + model.M₂ + model.M₃

export LangmuirS3