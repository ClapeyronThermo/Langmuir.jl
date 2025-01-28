"""
`LangmuirFreundlichS3(M₁, K₀₁, E₁, f₀₁, β₁, M₂, K₀₂, E₂, f₀₂, β₂, M₃, K₀₃, E₃, f₀₃, β₃)`

LangmuirFreundlichS3 <: IsothermModel

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
- `M₃::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₃::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₃::T`: Adsorption energy, `[J⋅mol⁻¹]`
- `f₀₃::T`: Surface heterogeneity parameter at T → ∞, `[-]`
- `β₃::T`: Surface heterogeneity coefficient, `[K]`

## Description

The Langmuir-Freundlich form of the isotherm for 3 sites is:

n = M₁ ⋅ K₁ ⋅ p₁^(f₀₁) ÷ (1 + K₁ ⋅ p₁^(f₀₁)) +  M₂ ⋅ K₂ ⋅ p₂^(f₀₂) ÷ (1 + K₂ ⋅ p₂^(f₀₂)) + M₃ ⋅ K₃ ⋅ p₃^(f₀₃) ÷ (1 + K₃ ⋅ p₃^(f₀₃))

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
@with_metadata struct LangmuirFreundlichS3{T} <: IsothermModel{T}
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
    (M₃::T, (0.0, Inf), "saturation loading 3")
    (K₀₃::T, (0.0, Inf), "Affinity parameter 3")
    (E₃::T, (-Inf, 0.0), "Energy parameter 3")
    (f₀₃::T, (0.0, Inf), "Surface heterogeneity parameter 3 at T → ∞")
    (β₃::T, (-Inf, Inf), "Surface heterogeneity coefficient 3")
end

function sp_res(model::LangmuirFreundlichS3, p, T)
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
    f₁ = model.f₀₁ - model.β₁/T 
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    f₂ = model.f₀₂ - model.β₂/T 
    K₃ = K₀₃*exp(-E₃/(Rgas(model)*T))
    f₃ = model.f₀₃ - model.β₃/T 
    return M₁*log1p(K₁*p^(f₁))/(f₁) + M₂*log1p(K₂*p^(f₂))/(f₂) + M₃*log1p(K₃*p^(f₃))/(f₃)
end

function loading(model::LangmuirFreundlichS3, p, T)
    M₁ = model.M₁
    K₀₁ = model.K₀₁
    E₁ = model.E₁
    M₂ = model.M₂
    K₀₂ = model.K₀₂
    E₂ = model.E₂
    M₃ = model.M₃
    K₀₃ = model.K₀₃
    E₃ = model.E₃
    f₁ = model.f₀₁ - model.β₁/T 
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    Kpf₁ = K₁*p^f₁
    f₂ = model.f₀₂ - model.β₂/T 
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    Kpf₂ = K₂*p^f₂
    f₃ = model.f₀₃ - model.β₃/T 
    K₃ = K₀₃*exp(-E₃/(Rgas(model)*T))
    Kpf₃ = K₃*p^f₃
    _1 = one(eltype(p))
    return M₁*Kpf₁/(_1 + Kpf₁) + M₂*Kpf₂/(_1 + Kpf₂) + M₃*Kpf₃/(_1 + Kpf₃)
end

saturated_loading(model::LangmuirFreundlichS3, T) = model.M₁ + model.M₂ + model.M₃

export LangmuirFreundlichS3