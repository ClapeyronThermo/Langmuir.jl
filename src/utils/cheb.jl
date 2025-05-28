## This is a vendored version of the code in EoSSuperancillaries.jl v1.4

struct ChebyshevRange{R,T}
    range::R
    coeffs::T
end

const ChebyshevRangeV64 = ChebyshevRange{Vector{Float64},Vector{Vector{Float64}}}
const ChebyshevRangeVec{T} = ChebyshevRange{Vector{T},Vector{Vector{T}}} where T

#overload of searchsortedfirst to also match the first element.
_searchsortedfirst(xx::Tuple,x) = _searchsortedfirst(SVector(xx),x)
_searchsortedfirst(xx,x) = Base.searchsortedfirst(xx,x) + isequal(x,first(xx))

function Base.show(io::IO,cheb::ChebyshevRange{R,T}) where {R,T}
    n = length(cheb.range) - 1
    min,max = extrema(cheb.range)
    print(io,typeof(cheb))
    print(io,"($n expansions from $min to $max)")
end

function Base.show(io::IO,::MIME"text/plain",cheb::ChebyshevRange{R,T}) where {R,T}
    println(io,typeof(cheb))
    n = length(cheb.range) - 1
    min,max = extrema(cheb.range)
    println(io,"  Expansions: $n")
    println(io,"  Minimum range: $min")
    print(io,"  Maximum range: $max")
end 

#evaluation of ranges of chebyshev coefficients
function cheb_eval(Base.@specialize(cheb), Base.@specialize(x̃))
    Cₙ = cheb.coeffs
    x̃range = cheb.range
    x̃min = first(x̃range)
    x̃max = last(x̃range)
    if !(x̃min <= x̃ <= x̃max)
        R = promote_type(eltype(x̃range), eltype(eltype(Cₙ)),typeof(x̃))
        #x is not in range
        return zero(R)/zero(R)
    end
    x̄,i = cheb_xrange(x̃range,x̃)
    Cₙi = Cₙ[i]
    return cheb_eval(Cₙi,x̄)
end

function cheb_xrange(x̃range,x̃)
        imax = _searchsortedfirst(x̃range, x̃)
        imin = imax - 1
        x̃minᵢ = x̃range[imin]
        x̃maxᵢ = x̃range[imax]
        i = imin
    
    #x̄ = (2*x̃ - (x̃maxᵢ + x̃minᵢ)) / (x̃maxᵢ - x̃minᵢ)
    x̄ = cheb_xtransform(x̃,x̃minᵢ,x̃maxᵢ)
    return x̄,i
end

cheb_xtransform(x,xmin,xmax) = (2*x - (xmax + xmin)) / (xmax - xmin)

function cheb_eval(cheb::AbstractVector{T},x::S) where {T,S}
    R = promote_type(T, S)
    l = length(cheb)
    l == 0 && return zero(R)
    l == 1 && return R(cheb[1])
    c0 = cheb[l - 1]
    c1 = cheb[l]
    for i in (l-2):-1:1
        c0, c1 = cheb[i] - c1, c0 + c1 * 2x
    end
    return R(c0 + c1 * x)
end

function ChebyshevRange(data::AbstractVector)
    first(data) isa AbstractDict || throw(ArgumentError("data must be a vector of dictionaries or named tuples"))
    n = length(data)
    c = Vector{Vector{Float64}}(undef,n)
    r = Vector{Float64}(undef,n + 1)
    r[1] = data[1][:xmin] #first value of the range, the minimum
    for i in 1:n
        data_i = data[i]
        coeff = data_i[:coef]
        xmax = data_i[:xmax]
        r[i+1] = xmax
        c[i] = coeff
    end
    return ChebyshevRange(r,c)
end

function chebnorm(M,coef::AbstractVector)
    coef_hi = @view coef[1:M]
    coef_lo = @view coef[(end - M + 1):end]
    return norm(coef_lo)/norm(coef_hi)
end

const CHEB_MAT_STORE = Dict{Tuple{Int,DataType},Matrix}()
const CHEB_POINT_STORE = Dict{Tuple{Int,DataType},Vector}()
const dCHEB_MAT_STORE = Dict{Tuple{Int,DataType},Matrix}()

function chebnorm(M,coef::Tuple)
   return chebnorm(M,SVector(coef))
end

chebpoints(N) = chebpoints(N,Float64)

function chebpoints(N,::Type{T}) where T
    n = T(N)
    return [cos(i*(pi/n)) for i in 0:N]
end

chebnorm(M) = Base.Fix1(chebnorm,M)

struct ChebTransformed{F,T}
    f::F
    xmin::T
    xmax::T
end

function (f̄::ChebTransformed{F,T})(x̄) where {F,T}
    x̄max,x̄min = f̄.xmax,f̄.xmin
    x = 0.5*x̄*(x̄max - x̄min) + 0.5*(x̄max + x̄min)
    convert(T,f̄.f(x))
end

function dct_mat(n::Int,type::Type{T}) where T
    L = zeros(T,n+1,n+1)
    return dct_mat!(L)
end

dct_mat(n::Int) = dct_mat(n,Float64)

function dct_mat(n::Val{N},type::Type{T}) where {N,T}
    L = @MMatrix zeros(T,N+1,N+1)
    return dct_mat!(L)
end

dct_mat(n::Val{N}) where {N} = dct_mat(n,Float64)

function dct_mat!(L::AbstractMatrix{T}) where T
    N = size(L,1) - 1
    N̄ = T(N)
    _2 = T(2)
    for j = 0:N
        for k = j:N
            p_j = (j == 0 || j == N) ? 2 : 1
            p_k = (k == 0 || k == N) ? 2 : 1
            L[j+1, k+1] = _2 / T(p_j * p_k * N) * cos(T(j * pi * k) / N̄)
            L[k+1, j+1] = L[j+1, k+1]
        end
    end
    return L
end

function cheb_interp(f,xmin::T,xmax::T,N::Int) where T
    f̄ = ChebTransformed(f,xmin,xmax)
    key = (N,typeof(xmin))
    x̄ = get!(CHEB_POINT_STORE,key) do
        chebpoints(N,T)::Vector{T}
    end
    M = get(CHEB_MAT_STORE,key) do
        dct_mat(N,T)::Matrix{T}
    end
    coef = M*map(f̄,x̄)
    return ChebyshevRange([xmin,xmax],[coef])
end

function merge_chebs!(left::ChebyshevRangeVec{T},right::ChebyshevRangeVec{T}) where T
    append!(left.range, @view(right.range[2:end]))
    left_coeffs = left.coeffs
    right_coeffs = right.coeffs
    for i in 1:length(right_coeffs)
        push!(left_coeffs, right_coeffs[i])
    end
    return left
end


"""
dyadic_splitting(f,N,xmin,xmax;M = 3,tol = 1.0e-12,max_refine_passes = 16,threaded = false)

Given a function `f`, with limits `xmin` and `xmax`,returns a `ChebishevRange`, with a N-range Chebishev interpolation at each interval.
The intervals are subdivided until `max_refine_passes` is reached or when the M-norm of chebishev coefficients is below `tol`.
"""
function dyadic_splitting(f,N,xmin,xmax;M = 3,tol = 1.0e-12,max_refine_passes = 16,threaded = false)
    #try to use 1 Chebyshev interpolation
    cheb = cheb_interp(f,xmin,xmax,N)
    if chebnorm(M,cheb.coeffs[1]) < tol
        return cheb
    end
    if max_refine_passes == 0
        @warn "max refine passes reached at x ∈ [$xmin,$xmax]"
        return cheb
    end
    #split domain
    xmid = 0.5*(xmin + xmax)
    left = dyadic_splitting(f,N,xmin,xmid,M = M,tol = tol,max_refine_passes = max_refine_passes - 1,threaded = threaded)
    right = dyadic_splitting(f,N,xmid,xmax,M = M,tol = tol,max_refine_passes = max_refine_passes - 1,threaded = threaded)
    return merge_chebs!(left,right)
end

#derivative of Chebyshev polynomial, modified from FastChebInterp.jl to just work in 1 dimension.
function cheb_deval(cheb::AbstractVector{T},x::S) where {T,S}
    R = promote_type(T, S)
    l = length(cheb)
    c₁ = cheb[1]
    l == 0 && return one(x) * zero(c₁), zero(R)/oneunit(x)
    l == 1 && return c₁ + one(x) * zero(c₁), zero(R)/oneunit(x)
    c₂ = cheb[2]
    l == 2 && return muladd(x, c₂, c₁), c₂*one(x)
    @inbounds cₙ₋₁ = cheb[1+(l-2)]
    @inbounds cₙ = cheb[1+(l-1)]
    bₖ′ = 2cₙ
    bₖ = muladd(x, bₖ′, cₙ₋₁)
    bₖ₊₁ = oftype(bₖ, cₙ)
    bₖ₊₁′ = zero(bₖ₊₁)
    for j = l-3:-1:1
        @inbounds cⱼ = cheb[1+j]
        bⱼ = muladd(2x, bₖ, cⱼ) - bₖ₊₁
        bⱼ′ = muladd(2x, bₖ′, 2bₖ) - bₖ₊₁′
        bₖ, bₖ₊₁ = bⱼ, bₖ
        bₖ′, bₖ₊₁′ = bⱼ′, bₖ′
    end
    return muladd(x, bₖ, c₁) - bₖ₊₁, muladd(x, bₖ′, bₖ) - bₖ₊₁′
end

function cheb_deval(Base.@specialize(cheb), Base.@specialize(x̃))
    Cₙ = cheb.coeffs
    x̃range = cheb.range
    x̃min = first(x̃range)
    x̃max = last(x̃range)
    if !(x̃min <= x̃ <= x̃max)
        R = promote_type(eltype(x̃range), eltype(eltype(Cₙ)),typeof(x̃))
        #x is not in range
        return zero(R)/zero(R), zero(R)/zero(R)
    end
    x̄,i = cheb_xrange(x̃range,x̃)
    Cₙi = Cₙ[i]
    return cheb_deval(Cₙi,x̄)
end

function derivate_cheb(ci::Vector{T}) where {T}
    cheb = Polynomials.ChebyshevT(ci)
    dcheb = Polynomials.derivative(cheb)
    return Polynomials.coeffs(dcheb)
end

function derivate_cheb(cheb::ChebyshevRangeVec{T}) where {T}
    dcoeffs = Vector{Vector{T}}(undef,length(cheb.coeffs))
    for i in 1:length(cheb.coeffs)
        dcoeffs[i] = derivate_cheb(cheb.coeffs[i])
    end
    return ChebyshevRange(cheb.range,dcoeffs)
end

#TODO: make this work
function integrate_cheb(ci::Vector{T},i0 = zero(T)) where {T}
    cheb = Polynomials.ChebyshevT{BigFloat}(ci)
    dcheb = Polynomials.integrate(cheb)
    return convert(Vector{T},Polynomials.coeffs(dcheb))
end

function integrate_cheb(cheb::ChebyshevRangeVec{T}) where {T}
    intcoeffs = Vector{Vector{T}}(undef,length(cheb.coeffs))
    i0 = zero(T)
    for i in 1:length(cheb.coeffs)
        icoeff = integrate_cheb(cheb.coeffs[i],i0)
        i0 = cheb_eval(icoeff,1.0) - i0
        
        intcoeffs[i] = icoeff
    end
    return ChebyshevRange(cheb.range,intcoeffs)
end
