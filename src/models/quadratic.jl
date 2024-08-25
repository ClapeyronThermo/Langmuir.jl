@with_metadata struct Quadratic{T} <: IsothermModel{T}
    (K₀a::T, (0.0, Inf), "Affinity parameter A")
    (K₀b::T, (0.0, Inf), "Affinity parameter B")
    (M::T, (0.0, Inf), "Saturation loading")
    (Ea::T, (-Inf, 0.0), "Energy parameter A")
    (Eb::T, (-Inf, 0.0), "Energy parameter B")
end

function sp_res(model::Quadratic, p, T)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    return M*log1p(p*(Ka + Kb*p))
end

function loading(model::Quadratic, p, T)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    return M*(Ka + 2*Kb*p)*p/(1 + p*(Ka + Kb*p))
end

function pressure_impl(model::Quadratic, Π, T, ::typeof(sp_res), approx)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    Kab = Ka/Kb
    return -0.5*Kab + sqrt(0.25*Kab*Kab + expm1(Π/M)/Kb)
end

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: Quadratic
    langmuir_model = x0_guess_fit(Langmuir,data)
    M,K₀,E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E
    #l*(1 + p*(Ka + Kb*p)) = M*(Ka*p + 2*Kb*p*p)
    #p*Ka*l + Kb*p*p*l - M*Ka*p - 2*M*Kb*p*p = -l
    #-p*Ka*l - Kb*p*p*l + M*Ka*p + 2*M*Kb*p*p = l
    #-Ka*(p*l) -Kb*(p*p*l) + M*Ka*(p) + 2*M*Kb*(p*p) = l

    l,p = data.l,data.p
    Kaneg,Kbneg,MKa,MKb2 = hcat(p .* l, p .* p .* l, p, p .*p ) \ l
    Ka = -Kaneg
    Kb = -Kbneg
    Ma = MKa/Ka
    Mb = 0.5*MKb2/Kb
    M = 0.5*(Ma + Mb)
    MK,K = hcat(p,-l .* p)\l
    M = MK/K
    _1 = one(eltype(langmuir_model))
    _0 = zero(eltype(langmuir_model))
    Quadratic(Ka, Kb, M, _0, _0)
end

export Quadratic
