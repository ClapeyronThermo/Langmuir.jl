"""
    `IsingS1(Mᵢ,Kᵢ₀,Eᵢ,Kₒ₀,Eₒ)`

    IsingS1 <: IsothermModel

`IsingS1(Mᵢ,Kᵢ₀,Eᵢ,Kₒ₀,Eₒ)` represents the single site Ising isotherm model.

## Inputs

- `Mᵢ::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `Kᵢ₀::T`: Affinity parameter I at T → ∞, `[Pa⁻¹]`
- `Eᵢ::T`: Adsorption energy I, `[J⋅mol⁻¹]`
- `Kₒ₀::T`: Affinity parameter O at T → ∞, `[Pa⁻¹]`
- `Eₒ::T`: Adsorption energy O, `[J⋅mol⁻¹]`


## Description

The IsingS1 equation is given by:

n = 

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀*exp(-E / (R * T))

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct IsingS1{T} <: IsothermModel{T}
    (Mᵢ ::T, (0.0, Inf), "saturation loading")
    (Kᵢ₀::T, (0.0, Inf), "affinity parameter I") #Using Inf cause trouble in bboxoptimize
    (Eᵢ::T, (-Inf, 0.0), "energy parameter I")
    (Kₒ₀::T, (0.0, Inf), "affinity parameter O") #Using Inf cause trouble in bboxoptimize
    (Eₒ::T, (-Inf, 0.0), "energy parameter O")
end

function loading(model::LangmuirS1, p, T)
    M = model.M
    Kₒ₀ = model.Kₒ₀
    Eₒ = model.Eₒ
    Kₒ = Kₒ₀*exp(-Eₒ/(Rgas(model)*T))
    Kᵢ₀ = model.Kᵢ₀
    Eᵢ = model.Eᵢ
    Kᵢ = Kᵢ₀*exp(-Eᵢ/(Rgas(model)*T))
    wᵢ = 0.5*(1-Kᵢ*p + ((1-Kᵢ*p)^2+4*Kₒ*p)^(0.5))

    return M * Kₒ *p / (wᵢ^2 + Kₒ*p)
end


