"""
`LangmuirFreundlichS4(M₁, K₀₁, E₁, f₀₁, β₁, M₂, K₀₂, E₂, f₀₂, β₂, M₃, K₀₃, E₃, f₀₃, β₃, M₄, K₀₄, E₄, f₀₄, β₄)`

LangmuirFreundlichS4 <: IsothermModel

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
- `M₄::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀₄::T`: Affinity parameter at T → ∞, `[1⋅Pa⁻¹]`
- `E₄::T`: Adsorption energy, `[J⋅mol⁻¹]`
- `f₀₄::T`: Surface heterogeneity parameter at T → ∞, `[-]`
- `β₄::T`: Surface heterogeneity coefficient, `[K]`

## Description

The Langmuir-Freundlich form of the isotherm for 4 sites is:

n = M₁ ⋅ K₁ ⋅ p₁^(f₀₁) ÷ (1 + K₁ ⋅ p₁^(f₀₁)) +  M₂ ⋅ K₂ ⋅ p₂^(f₀₂) ÷ (1 + K₂ ⋅ p₂^(f₀₂)) + M₃ ⋅ K₃ ⋅ p₃^(f₀₃) ÷ (1 + K₃ ⋅ p₃^(f₀₃)) + M₄ ⋅ K₄ ⋅ p₄^(f₀₄) ÷ (1 + K₄ ⋅ p₄^(f₀₄))

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
@with_metadata struct LangmuirFreundlichS4{T} <: IsothermModel{T}
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
    (M₄::T, (0.0, Inf), "saturation loading 4")
    (K₀₄::T, (0.0, Inf), "Affinity parameter 4")
    (E₄::T, (-Inf, 0.0), "Energy parameter 4")
    (f₀₄::T, (0.0, Inf), "Surface heterogeneity parameter 4 at T → ∞")
    (β₄::T, (-Inf, Inf), "Surface heterogeneity coefficient 4")
end


function sp_res(model::LangmuirFreundlichS4, p, T)
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
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    f₁ = model.f₀₁ - model.β₁/T 
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    f₂ = model.f₀₂ - model.β₂/T 
    K₃ = K₀₃*exp(-E₃/(Rgas(model)*T))
    f₃ = model.f₀₃ - model.β₃/T 
    K₄ = K₀₄*exp(-E₄/(Rgas(model)*T))
    f₄ = model.f₀₄ - model.β₄/T 
    return M₁*log1p(K₁*p^(f₁))/(f₁) + M₂*log1p(K₂*p^(f₂))/(f₂) + M₃*log1p(K₃*p^(f₃))/(f₃) + M₄*log1p(K₄*p^(f₄))/(f₄)
end

function loading(model::LangmuirFreundlichS4, p, T)
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
    f₁ = model.f₀₁ - model.β₁/T 
    K₁ = K₀₁*exp(-E₁/(Rgas(model)*T))
    Kpf₁ = K₁*p^f₁
    f₂ = model.f₀₂ - model.β₂/T 
    K₂ = K₀₂*exp(-E₂/(Rgas(model)*T))
    Kpf₂ = K₂*p^f₂
    f₃ = model.f₀₃ - model.β₃/T 
    K₃ = K₀₃*exp(-E₃/(Rgas(model)*T))
    Kpf₃ = K₃*p^f₃
    f₄ = model.f₀₄ - model.β₄/T 
    K₄ = K₀₄*exp(-E₄/(Rgas(model)*T))
    Kpf₄ = K₄*p^f₄
    _1 = one(eltype(p))
    return M₁*Kpf₁/(_1 + Kpf₁) + M₂*Kpf₂/(_1 + Kpf₂) + M₃*Kpf₃/(_1 + Kpf₃) + M₄*Kpf₄/(_1 + Kpf₄)
end

saturated_loading(model::LangmuirFreundlichS4, T) = model.M₁ + model.M₂ + model.M₃ + model.M₄

export LangmuirFreundlichS4
