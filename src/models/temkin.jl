struct TemkinApprox{T} <: IsothermModel{T}
    M::T
    K₀::T
    theta::T
    E::T
end

function sp_res(model::TemkinApprox, p, T)
    M, K₀, θ, E = model.M, model.K₀, model.theta, model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    Kp = K*p
    #is this abslopg1p alright?, removes errors when calling loading, but i don't know if its correct.
    return M*(log1p(abs(Kp)) + θ*(2*Kp + 1)/(2*(Kp + 1)*(Kp + 1)))
end

function loading(model::TemkinApprox, p, T)
    M, K₀, θ, E = model.M, model.K₀, model.theta, model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    Kp = K*p
    lm = Kp/(1 + Kp) #langmuir
    return M*lm*(1 + θ*lm*(lm - 1))
end

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: TemkinApprox
    #unilan ≈ langmuir with f = 1
    langmuir_model = x0_guess_fit(LangmuirS1,data)
    M,K₀,E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E
    return TemkinApprox(M, K₀,zero(eltype(langmuir_model)),E)
end

export TemkinApprox