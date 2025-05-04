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

include("ad.jl")
include("cheb.jl")
