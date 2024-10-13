#general multisite langmuir.

"""
    `Unilan(M, K₀, E)`

    Unilan <: IsothermModel

`Unilan(M, K₀, E, f)` represents the Unilan (UNIform distribution LANgmuir) isotherm model, which describes the adsorption of a gas on a solid surface.

## Inputs

- `M`::T: maximum loading capacity of the adsorbent, `[mol/kg]`
- `K₀`::T: equilibrium constant at zero coverage, `[1/Pa]`
- `E`::T: adsorption energy, `[J/mol]`
- `f`::T: heterogeneity of the adsorbent (no units)

## Description

The UNILAN equation is given by:

n = M * log((1 + K₀* exp(f) * p)/(1 + K₀ * exp(-f) * p)) / (2 * f)

where:
- n is the loading of the adsorbate on the adsorbent,
- M is the maximum loading capacity of the adsorbent,
- K₀ is the equilibrium constant at zero coverage,
- p is the pressure of the gas.
- f is the heterogeneity of the adsorbent. at the limit f -> 0, the langmuir isotherm is recovered.

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K₀ = exp(-E / (R * T))

where:
- R is the gas constant,
- T is the temperature.

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
