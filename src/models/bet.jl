struct BrunauerEmmettTeller{T} <: IsothermModel{T}
    K₀a::T
    K₀b::T
    M::T
    E::T
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

export BET
