
"""
    Toth <: IsothermModel

    Toth(M, K₀, f₀, β)

## Inputs

 - `M`::T: maximum loading capacity of the adsorbent, `[mol/kg]`
 - `K₀::T: equilibrium constant at zero coverage, `[1/Pa]`
 - `E`::T: adsorption energy, `[J/mol]`

## Description
K = K₀*exp(-E/(RT))
f = f₀ + β/T
nᵢ = M*K*P/(1 + (K*P)ᶠ)¹/ᶠ

"""
struct Toth{T} <: IsothermModel{T}
    M::T
    K₀::T
    f₀::T
    β::T
end


function sp_res(model::Toth, p, T)
#TODO (No analyical solution exist)

end


henry_coefficient(model::Toth, T) = #TODO
saturated_loading(model::Toth, T) = model.M #Some depend on T, some don't
sp_res_pressure_impl(model::Toth, Π, T) = #TODO


export Toth, henry_coefficient