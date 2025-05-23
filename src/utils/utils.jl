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
function unsafe_LU!(A::AbstractMatrix{T}, ipiv = Vector{BlasInt}(undef, min(size(A)...))) where {T}
    #check && LAPACK.chkfinite(A)
    # Extract values
    m, n = size(A)
    minmn = min(m,n)

    # Initialize variables
    info = 0
    #ipiv = Vector{BlasInt}(undef, minmn)
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if k < m
                amax = abs(A[k, k])
                for i = k+1:m
                    absi = abs(A[i,k])
                    if absi > amax
                        kp = i
                        amax = absi
                    end
                end
            end
            ipiv[k] = kp
            if !iszero(A[kp,k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k,i]
                        A[k,i] = A[kp,i]
                        A[kp,i] = tmp
                    end
                end
                # Scale first column
                Akkinv = inv(A[k,k])
                for i = k+1:m
                    A[i,k] *= Akkinv
                end
            elseif info == 0
                info = k
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end
        end
    end
    #check && checknonsingular(info, pivot)
    return LU(A, ipiv, convert(LinearAlgebra.BlasInt, info))
end

function unsafe_LU!(F::LU)
    return unsafe_LU!(F.factors,F.ipiv)
end

include("ad.jl")
include("cheb.jl")
