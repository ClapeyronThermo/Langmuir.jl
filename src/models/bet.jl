"""
    `BrunauerEmmettTeller(K₀a, K₀b, M, E)`

    BrunauerEmmettTeller <: IsothermModel

The Brunauer-Emmett-Teller (BET) isotherm model describes multilayer adsorption on a homogeneous surface. It is widely used to characterize adsorption behavior in mesoporous materials.

## Inputs

- `K₀a::T`: Affinity parameter A at T → ∞, `[1/Pa]`
- `K₀b::T`: Affinity parameter A at T → ∞, `[1/Pa]`
- `M::T`: Saturation loading, `[mol/kg]`
- `E::T`: Adsorption energy, `[J/mol]`

## Description

The BET equation is given by:

n = M × K₀a × p / (1 - K₀b × p) × (1 + (K₀a - K₀b) × p)


### Temperature dependence:
The affinity parameters `K₀a` and `K₀b` depend on temperature and are given by the following expressions:

K₀a = K₀a × exp(-E / (RT))

K₀b = K₀b × exp(-E / (RT))

Where:
- `R` is the gas constant,
- `T` is the absolute temperature.
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
