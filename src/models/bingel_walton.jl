"""
    `BingelWalton(M, K₀, A, E)`

    BingelWalton <: IsothermModel

`BingelWalton(M, K₀, A, E)` represents the Bingel & Walton isotherm model.

## Inputs

- `M::T`: Saturation loading, `[mol⋅kg⁻¹ ]`
- `K₀::T`: Intrinsic adsorption aﬃnity, `[1⋅Pa⁻¹]`
- `A::T`:  Clustering coeﬃcient describing strong adsorbate-adsorbate interactions, `[1⋅Pa⁻¹]`
- `E::T`: Adsorption energy, `[J⋅mol⁻¹]`

## Description

The Bingel & Walton equation is given by:

n = M⋅(1-exp(-(K+A)⋅p))÷(1+(A/K)⋅exp(-(K+A)⋅p))

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀⋅exp(-E÷(R⋅T))

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct BingelWalton{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading") 
    (K₀::T, (0.0, Inf), "affinity parameter") 
    (A::T, (0.0, Inf), "clustering coefficient") 
    (E::T, (-Inf, 0.0), "energy parameter")
end

#No analytical expression for the reduced grand potential 

function loading(model::BingelWalton, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    A = m0del.A
    K = K₀*exp(-E/(Rgas(model)*T))
    _1 = one(eltype(model))
    return M * (_1 - exp(-(K + A) * p)) / (_1 + (A/K) * exp(-(K + A) * p))
end

#optimizations for BingelWalton, not necessary, but improve performance
henry_coefficient(model::BingelWalton, T) = model.K₀*exp(-E/(Rgas(model)*T)) * model.M
saturated_loading(model::BingelWalton, T) = model.M #Some depend on T, some don't