"""
    `LangmuirS1(M, K₀, E)`

    LangmuirS1 <: IsothermModel

`LangmuirS1(M, K₀, E)` represents the single site Langmuir isotherm model.

## Inputs

- `M::T`: Saturation loading, `[mol⋅kg⁻¹]`
- `K₀::T`: Affinity parameter at T → ∞, `[Pa⁻¹]`
- `E::T`: Adsorption energy, `[J⋅mol⁻¹]`

## Description

The LangmuirS1 equation is given by:

n = (M * K * p) / (1 + K * p)

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀*exp(-E / (R * T))

Where:
- `R` is the universal gas constant, `[J⋅mol⁻¹⋅K⁻¹]`,
- `T` is the temperature, `[K]`.
"""

@with_metadata struct LangmuirS1{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "affinity parameter") #Using Inf cause trouble in bboxoptimize
    (E::T, (-Inf, 0.0), "energy parameter")
end

function sp_res(model::LangmuirS1, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    return M*log1p(K*p)
end


function loading(model::LangmuirS1, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))

    return M * K *p / (_1 + K*p)
end

#optimizations for LangmuirS1, not necessary, but improve performance
henry_coefficient(model::LangmuirS1, T) = model.M*model.K₀*exp(-model.E/(Rgas(model)*T))
saturated_loading(model::LangmuirS1, T) = model.M #Some depend on T, some don't
pressure_impl(model::LangmuirS1, Π, T,::typeof(sp_res), approx) = expm1(Π/model.M)/(model.K₀*exp(-model.E/(Rgas(model)*T)))

#TODO: include effects of temperature. at the moment, the fit procedure ignores temperature dependence.
#probably requires separating the models by temperature and linearizing K to obtain T-dependence.

function x0_guess_fit(::Type{T}, data::AdsIsoTData) where T <: LangmuirS1
    # use first two data points to get the slope

    #l = M*k*p/(1 + k*p)
    #l*(1 + k*p) = M*k*p
    #l + l*k*p = M*k*p
    #M*k*p - l*k*p = l
    #p*(Mk) - l*p(k) = l
    
    # Split data by temperature
    Ts, l_p = split_data_by_temperature(data)

    # Initialize vectors to store MK and K values
    MKs = Vector{eltype(Ts)}(undef, length(l_p))
    Ks = Vector{eltype(Ts)}(undef, length(l_p))

    # Perform the fitting for each (l, p) tuple
    for i in eachindex(l_p)
        l_min, p_min = l_p[i]
        MK, K = hcat(p_min, -l_min .* p_min) \ l_min
        MKs[i] = MK
        Ks[i] = K
    end

    M = sum(MKs./Ks)/length(Ks) #Mean of all values

    # log(K) = log(K0) - E/RT

    _1 = one(eltype(Ts))
    _1s = ones(eltype(Ts), length(Ts))

    if length(l_p) > 1
        logK, E = hcat(_1s, _1./ (Rgas(T).*Ts)) \ log.(Ks)
        K = exp(logK)
    else
        K = first(Ks)
        E = _1
    end
    
    return LangmuirS1(M, K, -E)
end

export LangmuirS1
