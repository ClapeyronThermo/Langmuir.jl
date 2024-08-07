Base.eltype(model::IsothermModel{T}) where T = T

function Base.eltype(::Type{M}) where M <: IsothermModel{T} where T
    return T
end

function (::Type{I})(p) where I <: IsothermModel
    return from_vec(I,p)
end

Rgas(model) = 8.31446261815324

"""
    model_length(model::IsothermModel)

Returns the amount of parameters necessary to instantiate a `model::IsothermModel`. For simple models, this is equivalent to `fieldcount(typeof(model))`, but composite models (like `MultiSite`), define their own implementation.
"""
#default.
model_length(::Type{T}) where T <: IsothermModel = _model_length(T)
model_length(model::IsothermModel) = model_length(typeof(model))

function _model_length(model::Type{T}) where T <: IsothermModel
    return fieldcount(T)
end

function from_vec(m::IsothermModel,x)
    return from_vec(typeof(m),x)
end

function from_vec(::Type{M},p::AbstractVector{K}) where {M <: IsothermModel,K}
    return M(ntuple(i -> p[i], model_length(M))...)
end

function from_vec(::Type{M},p::NTuple{N,K}) where {M <: IsothermModel,N,K}
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

#default: fill with ones
function x0_guess_fit(::Type{T}, data) where T <: IsothermModel
    p = pressure(data)
    _1 = one(eltype(p))
    _ones = ntuple(Returns(_1),model_length(T))
    _lb_bounds = 2 .* lower_bound(typeof(_1),T)
    v = max.(_ones,_lb_bounds)
    return from_vec(T,v)
end

"""
    isotherm_lower_bound(model::IsothermModel)
    isotherm_lower_bound(T,model::IsothermModel)
    isotherm_lower_bound(T,::Type{M}) where M <:IsothermModel

Returns the lower bound for the parameters of the isotherm model `model` of type `M`. with number type `T`, as a `Ntuple{model_length(M),T}`.
The default assumes that all parameters are nonnegative.
"""
function isotherm_lower_bound(model::IsothermModel)
    return isotherm_lower_bound(eltype(model),typeof(model))
end

function isotherm_lower_bound(::Type{T},m::IsothermModel) where T
    return isotherm_lower_bound(T,typeof(m))
end

function isotherm_lower_bound(::Type{T},::Type{M}) where T where M <: IsothermModel
    ntuple(Returns(zero(T)),model_length(M))
end

"""
    isotherm_upper_bound(model::IsothermModel)
    isotherm_upper_bound(T,model::IsothermModel)
    isotherm_upper_bound(T,::Type{M}) where M <:IsothermModel

Returns the upper bound for the parameters of the isotherm model `model` of type `M`. with number type `T`, as a `Ntuple{model_length(M),T}`.
The default assumes no upper bound for the parameters.
"""
function isotherm_upper_bound(model::T) where T <: IsothermModel
    return isotherm_upper_bound(eltype(model),typeof(model))
end

function isotherm_upper_bound(::Type{T},m::IsothermModel) where T
    return isotherm_upper_bound(T,typeof(m))
end

function isotherm_upper_bound(::Type{T},::Type{M}) where T where M <: IsothermModel
    ntuple(Returns(T(Inf)),model_length(M))
end

export isotherm_lower_bound, isotherm_upper_bound, model_length

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
include("Toth.jl")
include("multisite.jl")
include("interpolation.jl")
