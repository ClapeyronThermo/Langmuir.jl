"""
    Henry(Kh, E)

    Henry <: IsothermModel

## Inputs

- `Kh::T`: Affinity parameter, or Henry’s constant, `[mol/kg/Pa]`
- `E::T`: Adsorption energy, `[J/mol]`

## Description

The Henry isotherm model describes the adsorption of gases on solid surfaces at low pressures, where the amount of adsorbate adsorbed is directly proportional to the pressure of the gas. This model is typically valid in the low-pressure limit of adsorption, where the adsorption sites are far from saturation and the interactions between adsorbed molecules are negligible.

The adsorption behavior follows Henry's law:

n = Kh * p

The energy parameter `E` is related to Henry's constant `Kh` by the equation:

K = Kh*exp(-E / (R * T))

Where:
- `R` is the universal gas constant, `[J/mol/K]`,
- `T` is the temperature, `[K]`.

"""
@with_metadata struct Henry{T} <: IsothermModel{T}
    (Kh::T, (0.0, Inf), "affinity parameter")
    (E::T, (-Inf, 0.0), "energy parameter")
end

function sp_res(model::Henry, p, T)
    Kh = model.Kh
    K = Kh*exp(-E/(Rgas(model)*T))
    return K*p
end

function henry_coefficient(model::Henry, T)
    Kh = model.Kh
    K = Kh*exp(-E/(Rgas(model)*T))
    return K
end

pressure_impl(model::Henry,Π,T,::typeof(sp_res),approx) = Π/model.Kh*exp(-E/(Rgas(model)*T))

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: Henry
    # use first two data points to get the slope

    Ts, l_p = split_data_by_temperature(data)

    Ks = Vector{eltype(Ts)}(undef, length(l_p))

    for i in eachindex(l_p)
        l_i, p_i = l_p[i]
        K = p_i \ l_i
        Ks[i] = K
    end

    if length(l_p) > 1
        logK, E = hcat(_1s, _1./ (Rgas(T).*Ts)) \ log.(Ks)
        K = exp(logK)
    else
        K = first(Ks)
        E = _1
    end

    return Henry(K, -E)
end

export Henry
