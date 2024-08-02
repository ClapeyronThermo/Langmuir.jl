struct MultiSite{T,𝕀} <: IsothermModel{T}
    isotherms::𝕀
end

#we suppose that the isotherms tuple was processed before this
function _multisite(::Type{T},isotherms::𝕀) where {T,𝕀}
    return MultiSite{T,𝕀}(isotherms)
end

_multisite(isotherms::I) where I = _multisite(eltype(first(isotherms)),isotherms)

function MultiSite(m_first::IsothermModel,m_rest::Vararg{I}) where I <:IsothermModel
    
    return _multisite((m_first,m_rest...))
end
function model_length(model::MultiSite{T,I}) where {T,I}
    return model_length(typeof(model))
end

function model_length(::Type{MultiSite{T,I}}) where {T,I}
    return _model_length_multi(I)
end

function _model_length_multi(I) where {I}
    if @generated
        types = fieldtypes(I)
        N = mapreduce(model_length,+,types)
        return :(N)
    else
        types = fieldtypes(I)
        N = mapreduce(model_length,+,types)
        return N
    end
end

function Base.zero(m::MultiSite{T,I}) where {T,I}
    isotherms = m.isotherms
    zero_isotherms = zero.(isotherms)
    return _multisite(T,zero_isotherms)
end

function Base.zero(::Type{MultiSite{T,I}}) where {T,I}
    zero_isotherms = zero.(fieldtypes(I))
    return _multisite(T,zero_isotherms)
end

function Base.iszero(m::MultiSite{T,I}) where {T,I}
    zeros = iszero.(m.isotherms)
    return all(zeros)
end

function to_vec!(model::MultiSite,x)
    isotherms = m.isotherms
    n_begin = 1
    for i in 1:length(isotherms)
        model_i = isotherms[i]
        ni = model_length(model_i)
        n_end = n_begin + ni - 1
        vec_i = @view(x[n_begin:n_end])
        to_vec!(model_i,vec_i)
        n_begin = n_end + 1
    end
    return x
end

function from_vec(::Type{M},x) where {M::MultiSite{T,I},T,I}
    isotherms = zero(M)
    n_begin = 1
    for i in 1:length(isotherms)
        model_i = isotherms[i]
        ni = model_length(model_i)
        n_end = n_begin + ni - 1
        vec_i = @view(x[n_begin:n_end])
        new_model_i = from_vec(model_i,vec_i)
        isotherms = Base.setindex(isotherms,new_model_i,i)
        n_begin = n_end + 1
    end
    return MultiSite{T,I}(isotherms)
end

function loading(model::MultiSite,p,T)
    result = zero(Base.promote_eltype(model,p,T))
    for model in model.isotherms
        if !iszero(model)
            result += loading(model,p,T)
        end
    end
    return result
end

function sp_res(model::MultiSite,p,T)
    result = zero(Base.promote_eltype(model,p,T))
    for model in model.isotherms
        if !iszero(model)
            result += sp_res(model,p,T)
        end
    end
    return result
end

function henry_coefficient(model::MultiSite,T)
    result = zero(Base.promote_eltype(model,T))
    for model in model.isotherms
        if !iszero(model)
            result += henry_coefficient(model,T)
        end
    end
    return result
end

function saturated_loading(model::MultiSite,T)
    result = zero(Base.promote_eltype(model,T))
    for model in model.isotherms
        if !iszero(model)
            result += saturated_loading(model,T)
        end
    end
    return result
end
