"""
    loading(model::BrunauerEmmettTeller, p, T)

Calculates the loading `n` of adsorbate on the adsorbent based on the BET isotherm model.

## Inputs

- `model`::BrunauerEmmettTeller: A BET model instance containing the parameters K₀a, K₀b, M, and E.
- `p`::T: Partial pressure of the gas, `[Pa]`.
- `T`::T: Temperature, `[K]`.

## Outputs

- `n`::T: The loading of the adsorbate on the adsorbent, `[mol/kg]`.

## Description

This function computes the amount of gas adsorbed (`n`) based on the BET (Brunauer-Emmett-Teller) isotherm model. The BET model extends the Langmuir adsorption model to account for multilayer adsorption, which is important for modeling physical adsorption processes.

The function uses the following steps:
1. The equilibrium constants `Ka` and `Kb` for the first and subsequent layers are computed using the affinity parameters `K₀a` and `K₀b`, as well as the energy parameter `E`:
   
   Ka = K₀a * exp(-E / (R * T))

   Kb = K₀b * exp(-E / (R * T))

   where `R` is the universal gas constant.

2. The loading `n` is then calculated using the equation:

   n = M * (Ka * p) / ((1 - Kb * p) * (1 - Kb + Ka * p))

   where:
   - `M` is the saturation loading (maximum loading capacity of the first layer),
   - `p` is the partial pressure of the gas.

The function accounts for the adsorption in multiple layers, with the affinity parameters determining the strength of adsorption for each layer.
"""
@with_metadata struct BrunauerEmmettTeller{T} <: IsothermModel{T}
    (K₀a::T,  (0.0, Inf), "affinity parameter")
    (K₀b::T, (0.0, Inf), "affinity parameter")
    (M::T,  (0.0, Inf), "saturation loading")
    (E::T, (-Inf, 0.0), "energy parameter")
end

const BET = BrunauerEmmettTeller

function sp_res(model::BrunauerEmmettTeller, p, T)
    K₀a, K₀b, M, E = model.K₀a, model.K₀b, model.M, model.E
    Ka = K₀a*exp(-E/(Rgas(model)*T))
    Kb = K₀b*exp(-E/(Rgas(model)*T))
    Kap = Ka*p
    Kbp = Kb*p
    return M*log((1.0 - Kbp + Kap)/(1.0 - Kbp))
end

function loading(model::BrunauerEmmettTeller, p, T)
    K₀a, K₀b, M, E = model.K₀a, model.K₀b, model.M, model.E
    Ka = K₀a*exp(-E/(Rgas(model)*T))
    Kb = K₀b*exp(-E/(Rgas(model)*T))
    Kap = Ka*p
    Kbp = Kb*p
    return M*Kap/((1 - Kbp)*(1 - Kb + Kap))
end

export BET
