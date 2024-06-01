struct BrunauerEmmettTeller{T} <: IsothermModel{T}
    Ka::T
    Kb::T
    M::T
end

const BET = BrunauerEmmettTeller

function sp_res(model::BrunauerEmmettTeller,p)
    Ka,Kb,M = model.Ka,model.Kb,model.M
    Kap = Ka*p
    return M*log1p(Kap/(1 - Kap))
end

export BET