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

to_newton(f::F) where F = Base.Fix1(to_newton,f)

function to_halley(f,x)
    f,df,d2f = f∂f∂2f(f,x)
    return f,f/df,df/d2f
end

to_halley(f::F) where F = Base.Fix1(to_halley,f)