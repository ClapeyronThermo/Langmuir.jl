"""
    `LangmuirFreundlichS2(M₁, K₀₁, E₁, f₀₁, β₁, M₂, K₀₂, E₂, f₀₂, β₂)`

    LangmuirFreundlichS2 <: IsothermModel

## Inputs

- `M₁::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₁::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₁::T`: Adsorption energy, `[J⋅mol⁻¹]`
- `f₀₁::T`: Surface heterogeneity parameter at T → ∞, `[-]`
- `β₁::T`: Surface heterogeneity coefficient, `[K]`
- `M₂::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₂::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₂::T`: Adsorption energy, `[J⋅mol⁻¹]`
- `f₀₂::T`: Surface heterogeneity parameter at T → ∞, `[-]`
- `β₂::T`: Surface heterogeneity coefficient, `[K]`

## Description

The Langmuir-Freundlich form of the isotherm for 2 sites is:

n = M₁ ⋅ K₁ ⋅ p₁^(f₀₁) ÷ (1 + K₁ ⋅ p₁^(f₀₁)) +  M₂ ⋅ K₂ ⋅ p₂^(f₀₂) ÷ (1 + K₂ ⋅ p₂^(f₀₂))

Where:
- `n` is the loading of the adsorbate on the adsorbent,

### Temperature dependence:
The affinity parameter `Kᵢ` is temperature-dependent and can be expressed as:

Kᵢ = K₀ᵢ ⋅ exp(-Eᵢ ÷ (R ⋅ T))

The surface heterogeneity parameter `fᵢ` is also temperature-dependent and can be expressed as:

fᵢ = f₀ᵢ - βᵢ ÷ T

Where:
- `R` is the gas constant,
- `T` is the absolute temperature.

"""
@with_metadata struct LangmuirFreundlichS2{T} <: IsothermModel{T}
    (M₁::T, (0.0, Inf), "saturation loading 1")
    (K₀₁::T, (0.0, Inf), "Affinity parameter 1")
    (E₁::T, (-Inf, 0.0), "Energy parameter 1")
    (f₀₁::T, (0.0, Inf), "Surface heterogeneity parameter 1 at T → ∞")
    (β₁::T, (-Inf, Inf), "Surface heterogeneity coefficient 1")
    (M₂::T, (0.0, Inf), "saturation loading 2")
    (K₀₂::T, (0.0, Inf), "Affinity parameter 2")
    (E₂::T, (-Inf, 0.0), "Energy parameter 2")
    (f₀₂::T, (0.0, Inf), "Surface heterogeneity parameter 2 at T → ∞")
    (β₂::T, (-Inf, Inf), "Surface heterogeneity coefficient 2")
end


function sp_res(model::LangmuirFreundlichS2, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    f₁ = model.f₀₁ - model.β₁/T 
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    f₂ = model.f₀₂ - model.β₂/T 
    return M₁*log1p(K₁*p^(f₁))/(f₁) + M₂*log1p(K₂*p^(f₂))/(f₂)
end

function loading(model::LangmuirFreundlichS2, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    f₁ = model.f₀₁ - model.β₁/T 
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    Kpf₁ = K₁*p^f₁
    f₂ = model.f₀₂ - model.β₂/T 
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    Kpf₂ = K₂*p^f₂
    _1 = one(eltype(p))
    return M₁*Kpf₁/(_1 + Kpf₁) + M₂*Kpf₂/(_1 + Kpf₂)
end

saturated_loading(model::LangmuirFreundlichS2, T) = model.M₁ + model.M₂

export LangmuirFreundlichS2