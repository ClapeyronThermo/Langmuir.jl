struct LangmuirFreundlich{T} <: IsothermModel{T}
    M::T
    K₀::T
    E::T
    f::T
end


function sp_res(model::LangmuirFreundlich, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    return M*log1p(K*p^f)/f
end

function loading(model::LangmuirFreundlich, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    f = model.f
    K = K₀*exp(-E/(Rgas(model)*T))
    Kpf = K*p^f
    return M*Kpf/(1 + Kpf)
end

#optimizations for LangmuirFreundlich, not necessary, but improve performance
saturated_loading(model::LangmuirFreundlich, T) = model.M #Some depend on T, some don't

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: LangmuirFreundlich
    langmuir_model = x0_guess_fit(Langmuir,data)
    M, K₀, E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E
    _1 = one(M)
    return T(M, K₀, E, _1)
end

export LangmuirFreundlich, DualSiteLangmuirFreundlich
