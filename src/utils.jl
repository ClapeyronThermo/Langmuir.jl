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


"""
    f∂f(f,x)

returns f and ∂f/∂x evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""

@inline function f∂f(f::F, x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    out = f(ForwardDiff.Dual{T,R,1}(x, ForwardDiff.Partials((oneunit(R),))))
    return ForwardDiff.value(out),  ForwardDiff.extract_derivative(T, out)
end

f∂f(f::F) where F = Base.Fix1(f∂f,f)

"""
    f∂f∂2f(f,x)

returns f,∂f/∂x,and ∂²f/∂²x and evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
@inline function f∂f∂2f(f::F,x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    out = ForwardDiff.Dual{T,R,1}(x, ForwardDiff.Partials((oneunit(R),)))
    _f,_df = f∂f(f,out)
    fx = ForwardDiff.value(_f)
    dfx = ForwardDiff.partials(_f).values[1]
    d2fx = ForwardDiff.partials(_df).values[1]
    return (fx,dfx,d2fx)
end

f∂f∂2f(f::F) where F = Base.Fix1(f∂f∂2f,f)

"""
    fgradf2(f,x1,x2)

returns f and ∇f(x),using `ForwardDiff.jl`
"""

function fgradf2(f::F,x1::R1,x2::R2) where{F,R1<:Real,R2<:Real}
    y1,y2 = promote(x1,x2)
    return fgradf2(f,y1,y2)
end

@inline function fgradf2(f::F,x1::R,x2::R) where{F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    _1 = oneunit(R)
    _0 = zero(R)
    dual1 = ForwardDiff.Dual{T,R,2}(x1, ForwardDiff.Partials((_1,_0)))
    dual2 = ForwardDiff.Dual{T,R,2}(x2, ForwardDiff.Partials((_0,_1)))
    out = f(dual1,dual2)
    ∂out = ForwardDiff.partials(out)
    return ForwardDiff.value(out),SVector(∂out.values)
end

function to_newton(f,x)
    f,df = f∂f(f,x)
    return f,f/df
end
@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

_eltype(x::AbstractVector{T}) where T = T
_eltype(x::NTuple{N,T}) where {N,T} = T
function _eltype(x::T) where T
    eltyp = eltype(x)
    if eltyp === Any
        return Float64
    else
        return eltyp
    end
end

@inline function gradient(f::F, x) where {F}
    return ForwardDiff.gradient(f,x)
end