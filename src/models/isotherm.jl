#struct used to store collections of isotherms.
#mockup,.at the moment, im storing those as a tuple.
struct Isotherm{T,N,I} <: IsothermModel{T}
    isotherms::I
end

function Isotherm(isotherms,active = ntuple(Returns(true),length(isotherms)))
    T = base.promote_eltype(isotherms...)
    return Isotherm{T,length(isotherms),typeof(isotherms)}(isotherms,active)
end

const Isotherm1{T,I} = Isotherm{T,1,I} where {T,I}

function sp_res(modelx::Isotherm, p, T)
    result = zero(Base.promote_eltype(modelx,p,T))
    for model in modelx.isotherms
        !iszero(model)
            result += sp_res(model, p, T)
        end
    end
    return result
end

function loading(modelx::Isotherm, p, T)
    result = zero(Base.promote_eltype(modelx,p,T))
    for model in modelx.isotherms
        !iszero(model)
            result += loading(model, p, T)
        end
    end
    return result
end

function henry_coefficient(model::Isotherm, T)
    active = model.active
    result = zero(Base.promote_eltype(model,T))
    isotherms = model.isotherms
    for i in eachindex(active)
        if active
            result += henry_coefficient(isotherms[i], T)
        end
    end
    return result
end

function saturated_loading(model::Isotherm)
    active = model.active
    result = zero(Base.promote_eltype(model))
    isotherms = model.isotherms
    for i in eachindex(active)
        if active
            result += saturated_loading(isotherms[i])
        end
    end
    return result
end

#optimizations for an Isotherm containing only one isotherm model.
function sp_res(model::Isotherm1, p, T)
    return sp_res(model.isotherms[1], p, T)
end

function loading(model::Isotherm1, p, T)
    return loading(model.isotherms[1], p, T)
end

function henry_coefficient(model::Isotherm1, T)
    return henry_coefficient(model.isotherms[1], T)
end

function saturated_loading(model::Isotherm1)
    return saturated_loading(model.isotherms[1])
end

function sp_res_pressure_impl(model::Isotherm1, Π, T)
    return sp_res_pressure_impl(model.isotherms[1], Π, T)
end
