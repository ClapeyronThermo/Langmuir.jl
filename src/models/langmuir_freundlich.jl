"""
LangmuirFreundlich(M, K₀, E, f₀, β)

    LangmuirFreundlich <: IsothermModel

The Langmuir-Freundlich isotherm model is used to represent adsorption on heterogeneous surfaces. It combines elements of both the Langmuir and Freundlich models.

## Inputs

- `M::T`: Saturation loading, `[mol/kg]`
- `K₀::T`: Affinity parameter, `[1/Pa]`
- `E::T`: Adsorption energy, `[J/mol]`
- `f₀::T`: Surface heterogeneity parameter at high temperature, `[-]`
- `β::T`: Surface heterogeneity coefficient, `[K]`

## Description

The Langmuir-Freundlich isotherm can be used to model adsorption processes with both heterogeneous surface interactions and a saturation limit. The general form of the isotherm is:

q = M × log(1 + (K × p^f)) / f

Where:
- `q` is the loading of the adsorbate on the adsorbent,
- `M` is the saturation loading,
- `K` is the affinity parameter (temperature dependent),
- `p` is the pressure of the gas,
- `f` is the surface heterogeneity parameter, which describes how the adsorption sites deviate from uniformity.

### Temperature dependence:
The affinity parameter `K` is temperature-dependent and can be expressed as:

K = K₀ × exp(-E / (RT))

The surface heterogeneity parameter `f` is also temperature-dependent and can be expressed as:

f = f₀ - β / T

Where:
- `R` is the gas constant,
- `T` is the absolute temperature.

"""
@with_metadata struct LangmuirFreundlich{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "Affinity parameter")
    (E::T, (-Inf, 0.0), "Energy parameter")
    (f₀::T, (0.0, Inf), "Surface heterogeneity parameter at T → ∞")
    (β::T, (0.0, Inf), "Surface heterogeneity coefficient")
end


function sp_res(model::LangmuirFreundlich, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    f = model.f₀ - model.β/T 
    return M*log1p(K*p^f)/f
end

function loading(model::LangmuirFreundlich, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    f = model.f₀ - model.β/T 
    K = K₀*exp(-E/(Rgas(model)*T))
    Kpf = K*p^f
    _1 = one(eltype(p))
    return M*Kpf/(_1 + Kpf)
end

#optimizations for LangmuirFreundlich, not necessary, but improve performance
saturated_loading(model::LangmuirFreundlich, T) = model.M #Some depend on T, some don't

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: LangmuirFreundlich
    langmuir_model = x0_guess_fit(LangmuirS1,data)
    M, K₀, E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E
    _0 = zero(M)
    _1 = one(M)
    return T(M, K₀, E, _1, _0)
end

export LangmuirFreundlich
