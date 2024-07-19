struct QuadraticIsotherm{T} <: IsothermModel{T}
    K₀a::T
    K₀b::T
    M::T
    Ea::T
    Eb::T
end

function sp_res(model::QuadraticIsotherm, p, T)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    return M*log1p(p*(Ka + Kb*p))
end

function sp_res_pressure_impl(model::QuadraticIsotherm, Π, T)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    Kab = Ka/Kb
    return -0.5*Kab + sqrt(0.25*Kab*Kab + expm1(Π/M)/Kb)
end

export QuadraticIsotherm
