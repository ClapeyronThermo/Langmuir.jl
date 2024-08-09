struct MultiSite{T,𝕀} <: IsothermModel{T}
    isotherms::𝕀
end

isotherm_types(odel::MultiSite) = isotherm_types(typeof(model))

const __MultiSite{I,T} = MultiSite{T,T}

function isotherm_types(::Type{M}) where M <: MultiSite
    return fieldtypes(only(fieldtypes(M)))
end

function isotherm_types(::Type{__MultiSite{I}}) where {I}
    return fieldtypes(I)
end
#we suppose that the isotherms tuple was processed before this
function _multisite(::Type{T},isotherms::𝕀) where {T,𝕀}
    return MultiSite{T,𝕀}(isotherms)
end

_multisite(isotherms::I) where I = _multisite(eltype(first(isotherms)),isotherms)

function MultiSite(m_first::IsothermModel,m_rest::Vararg{I}) where I <:IsothermModel
    return _multisite((m_first,m_rest...))
end

function Base.getindex(m::MultiSite{T,𝕀}) where {T,𝕀}
    return m.isotherms[i]
end

"""

    MultiSite(m1::IsothermModel,isotherms::Vararg{IsothermModel}...)

given a list of isotherms, create a multisite isotherm model.

```julia
model1 = Langmuir(3.0,1.0,0.0)
model2 = Langmuir(3.0,0.9,3000.0)
double_site_langmuir = MultiSite(model1,model2) #create a multisite model with two langmuir isotherms

@assert loading(model1,1.0,3000.0) + loading(model2,1.0,3000.0) = loading(double_site_langmuir,1.0,3000.0)
```
"""

function model_length(::Type{MultiSite{T,I}}) where {T,I}
    return _model_length_multi(I)
end

function _model_length_multi(I::𝕀) where {𝕀}
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
    isotherms = model.isotherms
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

const _MultiSite{I,T} = MultiSite{T,I}

function from_vec(::Type{M},x,check) where M <: MultiSite
    MT = eltype(M)
    if MT === Any
        TT = _eltype(x)
    else
        TT = MT
    end
    iso_types = isotherm_types(M)
    isotherms = isotherm_zero.(TT,iso_types)
    n_begin = 1
    for i in 1:length(isotherms)
        model_i = isotherms[i]
        ni = model_length(model_i)
        n_end = n_begin + ni - 1
        vec_i = @view(x[n_begin:n_end])
        new_model_i = from_vec(model_i,vec_i,check)
        isotherms = Base.setindex(isotherms,new_model_i,i)
        n_begin = n_end + 1
    end
    return _multisite(isotherms)
end

function loading(model::MultiSite,p,T)
    result = zero(Base.promote_eltype(model,p,T))
    for model_i in model.isotherms
        if !iszero(model_i)
            result += loading(model_i,p,T)
        end
    end
    return result
end

function sp_res(model::M,p,T) where M <: MultiSite
    result = zero(Base.promote_eltype(model,p,T))
    for model_i in model.isotherms
        if !iszero(model_i)
            result += sp_res(model_i,p,T)
        end
    end
    return result
end

function henry_coefficient(model::MultiSite,T)
    result = zero(Base.promote_eltype(model,T))
    for model_i in model.isotherms
        if !iszero(model_i)
            result += henry_coefficient(model_i,T)
        end
    end
    return result
end

function saturated_loading(model::MultiSite,T)
    result = zero(Base.promote_eltype(model,T))
    for model_i in model.isotherms
        if !iszero(model)
            result += saturated_loading(model_i,T)
        end
    end
    return result
end

function isotherm_lower_bound(::Type{TT},::Type{MultiSite{T,I}}) where {TT,T,I}
    isotherms = fieldtypes(I)
    return tuplejoin(isotherm_lower_bound.(TT,isotherms))
end

function isotherm_upper_bound(::Type{TT},::Type{MultiSite{T,I}}) where {TT,T,I}
    isotherms = fieldtypes(I)
    return tuplejoin(isotherm_upper_bound.(TT,isotherms))
end

"""
    @MultiSite(isotherms)

Utility macro to build `MultiSite` types.

## Example:

```julia-repl
julia> v1 = @MultiSite{Langmuir,Langmuir} #abstract eltype
MultiSite{T, Tuple{Langmuir{T}, Langmuir{T}}} where T

julia> AdsorbedSolutionTheory.from_vec(v,[1,2,0,3,4,0])
MultiSite{Int64, Tuple{Langmuir{Int64}, Langmuir{Int64}}}((Langmuir{Int64}(1, 2, 0), Langmuir{Int64}(3, 4, 0)))

julia> v2 = @MultiSite{Langmuir,Langmuir}{Float64} #concrete eltype
MultiSite{Float64, Tuple{Langmuir{Float64}, Langmuir{Float64}}}

julia> AdsorbedSolutionTheory.from_vec(v2,[1,2,0,3,4,0])
MultiSite{Float64, Tuple{Langmuir{Float64}, Langmuir{Float64}}}((Langmuir{Float64}(1.0, 2.0, 0.0), Langmuir{Float64}(3.0, 4.0, 0.0)))
```
"""
macro MultiSite(isotherms)
    isotherms.head == :braces || throw(ArgumentError("invalid argument for `@MultiSite`: $isotherms"))
    IT = map(x-> Expr(:curly,x,:T),isotherms.args)
    I = Expr(:tuple,IT...)
    return quote
        MultiSite{T,Tuple{($I)...}} where T
    end |> esc
end

export MultiSite, @MultiSite
