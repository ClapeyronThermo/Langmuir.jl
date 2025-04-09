#general multisite langmuir.

"""
    `Unilan(M, K₀, E, f₀, β)`

    Unilan <: IsothermModel

## Inputs

- `M::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀::T`: Affinity parameter at T → ∞, [Pa⁻¹]`
- `E::T`: Adsorption energy, `[J⋅mol⁻¹]`
- `f₀::T`: Surface heterogeneity parameter at T → ∞, `[-]`
- `β::T`: Surface heterogeneity coefficient, `[K]`

## Description

The UNILAN equation is given by:

n = M * log((1 + K* exp(f) * p)/(1 + K * exp(-f) * p)) / (2 * f)

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀ × exp(-E / (R * T))

The surface heterogeneity parameter `f` is also temperature-dependent and can be expressed as:

f = f₀ - β / T

Where:
- `R` is the universal gas constant, `[J⋅mol⁻¹⋅K⁻¹]`,
- `T` is the temperature, `[K]`.
"""

@with_metadata struct Unilan{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "affinity parameter")
    (E::T, (-Inf, 0.0), "energy parameter")
    (f₀::T, (0.0, Inf), "surface heterogeneity parameter at T → ∞")
    (β::T, (0.0, Inf), "surface heterogeneity coefficient")
end

#b = K
function loading(model::Unilan, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    f = model.f₀ - model.β/T
    K = K₀*exp(-E/(Rgas(model)*T))
    Kfp1 = log1p(K*exp(f)*p)
    Kfp2 = log1p(K*exp(-f)*p)
    return M*(Kfp1 - Kfp2)/(2*f)
end

function sp_res(model::Unilan, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    f = model.f₀ - model.β/T
    K = K₀*exp(-E/(Rgas(model)*T))
    Li₂Kfp1 = PolyLog.reli2(-K*exp(-f)*p)
    Li₂Kfp2 = PolyLog.reli2(-K*exp(f)*p)
    (M/(2*f))*(Li₂Kfp1 - Li₂Kfp2)
end

#optimizations for unilan, not necessary, but improve performance
function henry_coefficient(model::Unilan, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    f = model.f₀ - model.β/T
    K = K₀*exp(-E/(Rgas(model)*T))
    (K*M/(2*f))*(exp(f) - exp(-f))
end

saturated_loading(model::Unilan, T) = model.M

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: Unilan
    #unilan ≈ langmuir with f = 1
    langmuir_model = x0_guess_fit(LangmuirS1,data)
    M, K₀, E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E
    return Unilan(M, K₀, E, one(eltype(langmuir_model)), zero(eltype(langmuir_model)))
end

export Unilan
