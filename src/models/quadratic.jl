struct QuadraticIsotherm{T} <: IsothermModel
    Ka::T
    Kb::T
    M::T
end

function sp_res(model::QuadraticIsotherm,p)
    Ka,Kb,M = model.Ka,model.Kb,model.M
    return M*log1p(p*(Ka+Kb*p))
end
