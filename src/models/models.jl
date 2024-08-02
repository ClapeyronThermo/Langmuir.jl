Base.eltype(model::IsothermModel{T}) where T = T

function Base.eltype(::Type{M}) where M <: IsothermModel{T} where T
    return T
end

function (::Type{I})(p) where I <: IsothermModel
    return from_vec(I,p)
end

Rgas(model) = 8.31446261815324

#default.
model_length(::Type{T}) where T <: IsothermModel = _model_length(T)
model_length(model::IsothermModel) = model_length(typeof(model))

function _model_length(model::Type{T}) where T <: IsothermModel
    return fieldcount(T)
end

function from_vec(m::IsothermModel,x) where K
    return from_vec(typeof(m),x)
end

function from_vec(::Type{M},Base.@specialize(p)) where M <: IsothermModel
    return M(ntuple(i -> p[i], model_length(M))...)
end

function to_vec!(model::IsothermModel,x)
    for i in 1:model_length(model)
        x[i] = getfield(model,i)
    end
    return x
end

function to_vec(model::IsothermModel)
    x = Vector{eltype(model)}(undef, model_length(model))
    to_vec!(model,x)
    return x
end

function to_tuple(model::IsothermModel)
    return ntuple(i -> getfield(model,i), model_length(model))
end

Base.zero(model::M) where M <: IsothermModel = Base.zero(M)

function Base.zero(model::Type{M}) where M <: IsothermModel{T} where T
    from_vec(M,ntuple(Returns(Base.zero(T)),model_length(model)))
end

#when T is not defined (zero(Langmuir))
function Base.zero(model::Type{M}) where M <: IsothermModel
    from_vec(M,ntuple(Returns(0.0),model_length(model)))
end

function Base.iszero(model::IsothermModel)
    result = true
    for i in 1:model_length(model)
        result &= iszero(getfield(model,i))
    end
    return result
end

function x0_guess_fit(::Type{T}, data) where T <: IsothermModel
    eltype = Base.promote_eltype(T, data)
end

export loading, henry_coefficient

include("freundlich.jl")
include("langmuir.jl")
include("langmuir_freundlich.jl")
include("redlich_peterson.jl")
include("sips.jl")
include("quadratic.jl")
include("bet.jl")
include("henry.jl")
include("temkin.jl")
include("unilan.jl")
include("toth.jl")
include("multisite.jl")
include("interpolation.jl")
