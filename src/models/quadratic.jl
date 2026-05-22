"""
    `Quadratic(K₀a, K₀b, M, Ea, Eb)`

    Quadratic <: IsothermModel

## Inputs

- `K₀a::T`: Affinity parameter A at T → ∞, `[1/Pa]`
- `K₀b::T`: Affinity parameter B at T → ∞, `[1/Pa^2]`
- `M::T`: Saturation loading, `[mol/kg]`
- `Ea::T`: Adsorption energy A, `[J/mol]`
- `Eb::T`: Adsorption energy B, `[J/mol]`

## Description

The `Quadratic` isotherm model is given by:

n = M × (Ka + 2Kb × p) × p / (1 + p × (Ka + Kb × p))

The model assumes that the affinity parameters `Ka` and `Kb` are temperature-dependent and follow the relation:

Ka = K₀a * exp(-Ea / (RT))

Kb = K₀b * exp(-Eb / (RT))

Where:
- `Ka` and `Kb` are the affinity parameters at temperature `T`,
- `R` is the gas constant,
- `T` is the absolute temperature.

"""
@with_metadata struct Quadratic{T} <: IsothermModel{T}
    (K₀a::T, (0.0, 1e2), "Affinity parameter A")
    (K₀b::T, (0.0, 1e2), "Affinity parameter B")
    (M::T, (0.0, 1e3), "Saturation loading")
    (Ea::T, (-2e5, 0.0), "Energy parameter A")
    (Eb::T, (-2e5, 0.0), "Energy parameter B")
end

function sp_res(model::Quadratic, p, T)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    return M*log1p(p*(Ka + Kb*p))
end

function loading(model::Quadratic, p, T)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    _1 = one(eltype(p))
    return M*(Ka + 2.0*Kb*p)*p/(_1 + p*(Ka + Kb*p))
end

henry_coefficient(model::Quadratic, T) = model.M*model.K₀a*exp(-model.Ea/(Rgas(model)*T))

function pressure_impl(model::Quadratic, Π, T, ::typeof(sp_res))
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    Kab = Ka/Kb
    return -0.5*Kab + sqrt(0.25*Kab*Kab + expm1(Π/M)/Kb)
end

saturated_loading(model::Quadratic, T) = 2*model.M

function x0_guess_fit(::Type{Q}, data::AdsIsoTData) where Q <: Quadratic

    # Simplified approach: Use Langmuir fit as initial guess for M,K0a,Ea
    langmuir_model = x0_guess_fit(LangmuirS1, data)
    M = langmuir_model.M*0.5
    K₀a = langmuir_model.K₀*2
    Ea = langmuir_model.E    
    p = pressure(data)
    T = temperature(data)
    l = loading(data)
    RT_inv = @. 1/(Rgas(1.0)*T)
    Ka = @. K₀a*exp(-Ea*RT_inv)

    #=
    l*(_1 + p*(Ka + Kb*p)) = M*(Ka + 2.0*Kb*p)*p
    l*(_1 + p*Ka + Kb*p*p) = M*(Ka*p * Kb*p*p)
    l + l*p*Ka + l*Kb*p*p - M*Ka*p - 2*M*Kb*p*p = 0
    l*Kb*p*p - 2*M*Kb*p*p = -(l + l*p*Ka - M*Ka*p)
    l*Kb - 2*M*Kb= -(l + l*p*Ka - M*Ka*p) ./ p*p
    Kb= -(l + l*p*Ka - M*Ka*p) ./ p*p / (l - 2M)
    =#
    Kb =  @. (-l/(p*p) - (l - M)*Ka/p) /(l - 2M)
    i = findall(x -> isfinite(x) && x > 0,Kb)
    Kb = Kb[i]
    RT_inv = RT_inv[i]
    Kb .= clamp.(Kb, 1e-30, 1e15)
    logKb = Kb
    logKb .= log.(Kb)
    _1s  = ones(eltype(RT_inv), length(RT_inv))
    coeffs = hcat(_1s, RT_inv) \ logKb
    logK₀b = coeffs[1]
    slope = coeffs[2]
    K₀b = exp(logK₀b)
     Eb = -slope  # Since slope = -E, we have E = -slope
    @show M,Ea,K₀a
   
    return Q(K₀a, K₀b, M, Ea, Eb)
end

export Quadratic
