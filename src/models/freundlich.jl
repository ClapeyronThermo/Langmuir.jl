struct Freundlich{T} <: IsothermModel{T}
    M::T
    f::T
end

function sp_res(model::Freundlich,p)
    M, f = model.M, model.f
    return M*p^f/f
end

function pressure_impl(model::Freundlich,Π,T,::typeof(sp_res),approx)
    M, f = model.M, model.f
    v = 1/f
    return (Π/(M*v))^v
end

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: Freundlich
    l,p = data.l,data.p
    #l = M*p^f
    #log(l) = log(M) + f*log(p)
    logp,logl = log.(p),log.(l)
    _1 = one.(p)
    logM,f = hcat(_1,logp) \ logl
    M = exp(logM)
    return T(M,f)
end

export Freundlich
