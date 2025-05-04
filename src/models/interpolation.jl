#struct to cache isotherm loading and sp_res
struct CachedIsotherm{T,I} <: IsothermModel{T}
    pmax::T
    temp::T
    isotherm::I
    sp_res_interpolator::ChebyshevRangeVec{T}
    loading_interpolator::ChebyshevRangeVec{T}
end

function CachedIsotherm(isotherm::I,pmax,T) where I
    TT = Base.promote_eltype(isotherm, pmax, T)
    if requires_integration_sp_res(isotherm) == Val{true}()
        f_loading(x) = loading(isotherm, x, T)
        loading_interpolator = dyadic_splitting(f_loading, 16, eps(TT), TT(pmax), max_refine_passes = 32)
        sp_res_interpolator = integrate_cheb(loading_interpolator)
        return CachedIsotherm{TT,I}(pmax,T,isotherm,sp_res_interpolator,loading_interpolator)
    else
        f_sp_res(x) = sp_res(isotherm, x, T)
        sp_res_interpolator = dyadic_splitting(f_sp_res, 16, eps(TT), TT(pmax),max_refine_passes = 32)
        loading_interpolator = derivate_cheb(sp_res_interpolator)
        return CachedIsotherm{TT,I}(pmax,T,isotherm,sp_res_interpolator,loading_interpolator)
    end
end

function Base.show(io::IO, ::MIME"text/plain", model::CachedIsotherm{T,I}) where {T,I}
    print(io,typeof(model))
    print(io,"(T = $(model.temp), pmax = $(model.pmax))") 
end

function Base.show(io::IO, model::CachedIsotherm{T,I}) where {T,I}
    print(io,typeof(model))
    print(io,"(T = $(model.temp), pmax = $(model.pmax))") 
end

function sp_res(model::CachedIsotherm, p, T)
    if p > model.pmax
        #fall back to the original model
        return sp_res(model.isotherm, p, T)
    else
        return cheb_eval(model.sp_res_interpolator, p)
    end
end

function loading(model::CachedIsotherm, p, T)
    pmax = model.pmax
    if p > pmax
        #fall back to the original model
        return loading(model.isotherm, p, T)
    else
        return cheb_eval(model.loading_interpolator, p)
    end
end

function henry_coefficient(model::CachedIsotherm, T)
    pmin = first(model.loading_interpolator.range)
    if pmin > 10eps(eltype(pmin))
        return henry_coefficient(model.isotherm, T)
    end
    cheb_zero_p = model.loading_interpolator.coeffs[1]
    _,dldp = cheb_deval(cheb_zero_p, -oneunit(T))
    return dldp
end
saturated_loading(model::CachedIsotherm, T) = saturated_loading(model.isotherm, T)

export CachedIsotherm