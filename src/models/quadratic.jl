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

function sp_res_pressure_impl(model::QuadraticIsotherm, Π, T, approx)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    Kab = Ka/Kb
    return -0.5*Kab + sqrt(0.25*Kab*Kab + expm1(Π/M)/Kb)
end

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: QuadraticIsotherm
    langmuir_model = x0_guess_fit(Langmuir,data)
    M,K₀,E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E
    _1 = one(eltype(langmuir_model))
    _0 = zero(eltype(langmuir_model))
    QuadraticIsotherm(M, K₀, _1, E, _0)
end

export QuadraticIsotherm
