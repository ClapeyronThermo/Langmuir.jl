struct QuadraticIsotherm{T} <: IsothermModel
    Ka::T
    Kb::T
    M::T
end

function sp_res(model::QuadraticIsotherm,p)
    Ka,Kb,M = model.Ka,model.Kb,model.M
    return M*log1p(p*(Ka+Kb*p))
end

function sp_res_pressure_impl(model::QuadraticIsotherm,p)
    Ka,Kb,M = model.Ka,model.Kb,model.M
    Kab = Ka/Kb
    return -0.5*Kab + sqrt(0.25*Kab*Kab + expm1(q/M)/Kb)
end

export QuadraticIsotherm
