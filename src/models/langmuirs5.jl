"""
    `LangmuirS5(M₁, K₀₁, E₁, M₂, K₀₂, E₂, M₃, K₀₃, E₃, M₄, K₀₄, E₄, M₅, K₀₅, E₅)`

    LangmuirS5 <: IsothermModel

`LangmuirS5(M₁, K₀₁, E₁, M₂, K₀₂, E₂, M₃, K₀₃, E₃, M₄, K₀₄, E₄, M₅, K₀₅, E₅)` represents the 5 site Langmuir isotherm model.

## Inputs

- `M₁::T`: Saturation loading 1, `[mol⋅kg⁻¹]`
- `K₀₁::T`: Affinity parameter 1 at T → ∞, `[1⋅Pa⁻¹]`
- `E₁::T`: Adsorption energy 1, `[J⋅mol⁻¹]`
- `M₂::T`: Saturation loading 2, `[mol⋅kg⁻¹]`
- `K₀₂::T`: Affinity parameter 2 at T → ∞, `[1⋅Pa⁻¹]`
- `E₂::T`: Adsorption energy 2, `[J⋅mol⁻¹]`
- `M₃::T`: Saturation loading 3, `[mol⋅kg⁻¹]`
- `K₀₃::T`: Affinity parameter 3 at T → ∞, `[1⋅Pa⁻¹]`
- `E₃::T`: Adsorption energy 3, `[J⋅mol⁻¹]`
- `M₄::T`: Saturation loading 4, `[mol⋅kg⁻¹]`
- `K₀₄::T`: Affinity parameter 4 at T → ∞, `[1⋅Pa⁻¹]`
- `E₄::T`: Adsorption energy 4, `[J⋅mol⁻¹]`
- `M₅::T`: Saturation loading 5, `[mol⋅kg⁻¹]`
- `K₀₅::T`: Affinity parameter 5 at T → ∞, `[1⋅Pa⁻¹]`
- `E₅::T`: Adsorption energy 5, `[J⋅mol⁻¹]`

## Description

The LangmuirS5 equation is given by:

n = (M₁ ⋅ K₁ ⋅ p) ÷ (1 + K₁ ⋅ p) + (M₂ ⋅ K₂ ⋅ p) ÷ (1 + K₂ ⋅ p) + (M₃ ⋅ K₃ ⋅ p) ÷ (1 + K₃ ⋅ p) + (M₄ ⋅ K₄ ⋅ p) ÷ (1 + K₄ ⋅ p) + (M₅ ⋅ K₅ ⋅ p) ÷ (1 + K₅ ⋅ p)

The adsorption energy Eᵢ is related to the equilibrium constant K₀ᵢ by the equation:      

Kᵢ = K₀ᵢ ⋅ exp(-Eᵢ ÷ (R ⋅ T))

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct LangmuirS5{T} <: IsothermModel{T}
    (M₁::T, (0.0, Inf), "saturation loading 1")
    (K₀₁::T, (0.0, Inf), "affinity parameter 1") 
    (E₁::T, (-Inf, 0.0), "energy parameter 1")
    (M₂::T, (0.0, Inf), "saturation loading 2")
    (K₀₂::T, (0.0, Inf), "affinity parameter 2") 
    (E₂::T, (-Inf, 0.0), "energy parameter 2")
    (M₃::T, (0.0, Inf), "saturation loading 3")
    (K₀₃::T, (0.0, Inf), "affinity parameter 3") 
    (E₃::T, (-Inf, 0.0), "energy parameter 3")
    (M₄::T, (0.0, Inf), "saturation loading 4")
    (K₀₄::T, (0.0, Inf), "affinity parameter 4") 
    (E₄::T, (-Inf, 0.0), "energy parameter 4")
    (M₅::T, (0.0, Inf), "saturation loading 5")
    (K₀₅::T, (0.0, Inf), "affinity parameter 5") 
    (E₅::T, (-Inf, 0.0), "energy parameter 5")
end

function sp_res(model::LangmuirS5, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    M₃ = model.M₃
    K₀₃ = model.K₀₃
    E₃ = model.E₃
    M₄ = model.M₄
    K₀₄ = model.K₀₄
    E₄ = model.E₄
    M₅ = model.M₅
    K₀₅ = model.K₀₅
    E₅ = model.E₅
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    K₃ = K₀₃*exp(-E₃/(Rgas(model)*T))
    K₄ = K₀₄*exp(-E₄/(Rgas(model)*T))
    K₅ = K₀₅*exp(-E₅/(Rgas(model)*T))
    return M₁*log1p(K₁*p)+M₂*log1p(K₂*p)+M₃*log1p(K₃*p)+M₄*log1p(K₄*p)+M₅*log1p(K₅*p)
end

function loading(model::LangmuirS5, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    M₃ = model.M₃
    K₀₃ = model.K₀₃
    E₃ = model.E₃
    M₄ = model.M₄
    K₀₄ = model.K₀₄
    E₄ = model.E₄
    M₅ = model.M₅
    K₀₅ = model.K₀₅
    E₅ = model.E₅
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    K₃ = K₀₃*exp(-E₃/(Rgas(model)*T))
    K₄ = K₀₄*exp(-E₄/(Rgas(model)*T))
    K₅ = K₀₅*exp(-E₅/(Rgas(model)*T))
    _1 = one(eltype(model))
    return M₁ * K₁ *p / (_1 + K₁*p)+ M₂ * K₂ *p / (_1 + K₂*p)+ M₃ * K₃ *p / (_1 + K₃*p) + M₄ * K₄ *p / (_1 + K₄*p) + M₅ * K₅ *p / (_1 + K₅*p)
end

henry_coefficient(model::LangmuirS5, T) = model.M₁*model.K₀₁*exp(-model.E₁/(Rgas(model)*T))+model.M₂*model.K₀₂*exp(-model.E₂/(Rgas(model)*T))+model.M₃*model.K₀₃*exp(-model.E₃/(Rgas(model)*T))+model.M₄*model.K₀₄*exp(-model.E₄/(Rgas(model)*T))+model.M₅*model.K₀₅*exp(-model.E₅/(Rgas(model)*T))
saturated_loading(model::LangmuirS5, T) = model.M₁ + model.M₂ + model.M₃ + model.M₄ + model.M₅

export LangmuirS5