#=
structs ported from Clapeyron.jl, and other utils
=#

struct FractionVector{T,V} <: AbstractVector{T}
    vec::V
    val::T
end
#=

Fraction Vector
useful when expressing fractions in n-1 components.
the last component is calculated at build time.
it allocates less than creating a new vector or appending.
=#
##
function FractionVector(v::AbstractVector)
    a = Base.one(eltype(v))
    # any(x->x<0,v) && throw(DomainError(v,"all elements of a fraction vector should be positive."))
    a -=Base.sum(v)
    # a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
    return FractionVector(v,a) 
end

function FractionVector(v::Real)
    a = Base.one(v) 
    # (v < zero(v)) && throw(DomainError(v,"all elements of a fraction vector should be positive."))
    a -= v
    # a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
    return FractionVector(v,a)
end


@inline Base.eltype(v::FractionVector{T}) where T = T
@inline Base.length(v::FractionVector)::Int = Base.length(v.vec) + 1

@inline function Base.length(v::FractionVector{T,<:Real})::Int where T
    return 2
end

@inline Base.size(v::FractionVector) = (length(v),)
@inline function Base.getindex(v::FractionVector,i::Int)
    @boundscheck checkbounds(v, i)
    if length(v)==i
         return v.val
    else
        return v.vec[i]
    end
end

@inline function Base.getindex(v::FractionVector{T,<:Real},i::Int) where T
    @boundscheck checkbounds(v, i)
    return ifelse(i==1,v.vec,v.val)
end

Base.IndexStyle(::Type{<:FractionVector}) = IndexLinear()
