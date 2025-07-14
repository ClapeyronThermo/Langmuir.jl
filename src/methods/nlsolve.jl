#nlsolve functionality
"""
    function nlsolve(f!,x0,method=TrustRegion(Newton(), Dogleg()), options=NEqOptions(),chunk = ForwardDiff.Chunk{2}())


Given a function `f!(result,x)` that returns a system of equations,
`nlsolve(f!,x0)` returns a `NLSolvers.ConvergenceInfo` struct that contains the results of the non-linear solving procedure.

Uses `NLSolvers.jl` as backend, the jacobian is calculated with `ForwardDiff.jl`, with the specified `chunk` size

To obtain the underlying solution vector, use [`x_sol`](@ref)

To see available solvers and options, check `NLSolvers.jl`
"""
function nlsolve(f!,x0,method=TrustRegion(Newton(), Dogleg()),options=NEqOptions(),chunk = ForwardDiff.Chunk{2}())
    vector_objective = autoVectorObjective(f!,x0,chunk)
    nl_problem = NEqProblem(vector_objective; inplace = _inplace(x0))
    return nlsolve(nl_problem, x0,method, options)
end

function nlsolve(nl_problem::NEqProblem,x0,method =TrustRegion(Newton(), NWI()),options=NEqOptions())
    return NLSolvers.solve(nl_problem, x0,method, options)
end

function autochunk(x)
    k = ForwardDiff.pickchunksize(length(x))
    return ForwardDiff.Chunk{k}()
end


function autoVectorObjective(f!,x0,chunk)
    Fcache = x0 .* false
    jconfig = ForwardDiff.JacobianConfig(f!,x0,x0,chunk)
    function j!(J,x)
        ForwardDiff.jacobian!(J,f!,Fcache,x,jconfig)
        J
    end
    function fj!(F,J,x)
        ForwardDiff.jacobian!(J,f!,F,x,jconfig)
        F,J
    end
    function jv!(x)
        return nothing
    end
    return NLSolvers.VectorObjective(f!,j!,fj!,jv!)
end

_inplace(x0) = true
_inplace(x0::SVector) = false
_inplace(x0::Number) = false
function autoVectorObjective(f!,x0::StaticArrays.SVector{2,T},chunk) where T
    f(x) = f!(nothing,x) #we assume that the F argument is unused in static arrays
    j(J,x) = ForwardDiff.jacobian(f,x)
    fj(F,J,x) = FJ_ad(f,x)
    return NLSolvers.VectorObjective(f!,j,fj,nothing)
end

function autoVectorObjective(f!,x0::StaticArrays.SVector{3,T},chunk) where T
    f(x) = f!(nothing,x) #we assume that the F argument is unused in static arrays
    j(J,x) = ForwardDiff.jacobian(f,x)
    fj(F,J,x) = FJ_ad(f,x)
    return NLSolvers.VectorObjective(f!,j,fj,nothing)
end

function autoVectorObjective(f!,x0::StaticArrays.SVector,chunk)
    f(x) = f!(nothing,x) #we assume that the F argument is unused in static arrays
    j(J,x) = ForwardDiff.jacobian(f,x)
    fj(F,J,x) = FJ_ad(f,x)
    return NLSolvers.VectorObjective(f,j,fj,nothing)
end

#= only_fj!: NLsolve.jl legacy form:

function only_fj!(F, J, x)
    # shared calculations begin
    # ...
    # shared calculation end
    if !(F == nothing)
        # mutating calculations specific to f! goes here
    end
    if !(J == nothing)
        # mutating calculations specific to j! goes
    end
end
=#
function only_fj!(fj!::T) where T
    function _f!(F,x)
        fj!(F,nothing,x)
        F
    end

    function _fj!(F,J,x)
        fj!(F,J,x)
        F,J
    end

    function _j!(J,x)
        fj!(nothing,J,x)
        J
    end
    return NLSolvers.VectorObjective(_f!,_j!,_fj!,nothing) |> NEqProblem
end

function ADScalarObjective(f,x0::AbstractArray,chunk = autochunk(x0))
    Hres = DiffResults.HessianResult(x0)
    function _g(df,x,Gresult)
        ForwardDiff.gradient!(Gresult,f,x)
        df .= DiffResults.gradient(Gresult)
        df
    end
    
    function _fg(df,x,Gresult)
        ForwardDiff.gradient!(Gresult,f,x)
        df .= DiffResults.gradient(Gresult)
        fx = DiffResults.value(Gresult)
        return fx,df
    end

    function _fgh(df,d2f,x,Hresult)
        ForwardDiff.hessian!(Hresult,f,x)
        d2f .= DiffResults.hessian(Hresult)
        df .= DiffResults.gradient(Hresult)
        fx = DiffResults.value(Hresult)
        return fx,df,d2f
    end

    function h(d2f,x)
        ForwardDiff.hessian!(d2f,f,x)
        d2f
    end
    g(df,x) = _g(df,x,Hres)
    fg(df,x) = _fg(df,x,Hres)
    fgh(df,d2f,x) = _fgh(df,d2f,x,Hres)
    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end


function ADScalarObjective(f,x0::Number,autochunk)
    function g(x)
        return derivative(f,x)
    end

    function fg(∂fx,x)
        ∂fx,x = f∂f(f,x)
        return ∂fx,x
    end

    function fgh(∂fx,∂2fx,x)
        fx,∂fx,∂2fx = f∂f∂2f(f,x)
        return fx,∂fx,∂2fx
    end

    function h(∂2fx,x)
        ∂2fx = derivative(g,x)
        return ∂2fx
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end
#uses brent, the same default that Optim.jl uses
function optimize(f,x0::NTuple{2,T},method=BrentMin(),options=OptimizationOptions()) where T<:Real
    scalarobj = ADScalarObjective(f,first(x0),nothing)
    optprob = NLSolvers.OptimizationProblem(obj = scalarobj,bounds = x0, inplace=false)
    res = NLSolvers.solve(optprob,method,options)
    return res
end
#general one, with support for ActiveBox
function optimize(f,x0,method=LineSearch(Newton()),options=OptimizationOptions();bounds = nothing)
    scalarobj = ADScalarObjective(f,x0,autochunk)
    optprob = NLSolvers.OptimizationProblem(obj = scalarobj,inplace = _inplace(x0),bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(optprob::NLSolvers.OptimizationProblem,x0,method=LineSearch(Newton()),options=OptimizationOptions();bounds = nothing)
    return NLSolvers.solve(optprob,x0,method,options)
end
#build scalar objective -> Optimization Problem
function optimize(scalarobj::ScalarObjective,x0,method=LineSearch(Newton()),options=OptimizationOptions();bounds = nothing)
    optprob = NLSolvers.OptimizationProblem(obj = scalarobj,inplace = _inplace(x0),bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(f,x0,method::NLSolvers.NelderMead,options=OptimizationOptions();bounds = nothing)
    scalarobj = ScalarObjective(f = f)
    optprob = NLSolvers.OptimizationProblem(obj = scalarobj,inplace = _inplace(x0),bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

x_minimum(res::NLSolvers.ConvergenceInfo) = res.info.minimum
#for BrentMin (should be fixed at NLSolvers 0.3)
x_minimum(res::Tuple{<:Number,<:Number}) = last(res)

struct RestrictedLineSearch{F,LS} <: NLSolvers.LineSearcher
    f::F #function that restricts the line search
    ls::LS #actual line search
end


function NLSolvers.find_steplength(mstyle::NLSolvers.MutateStyle, ls::RestrictedLineSearch{F,LS}, φ::T, λ) where {F,LS,T}
    λr = ls.f(φ,λ)
    NLSolvers.find_steplength(mstyle, ls.ls, φ, λ)
end

function ls_restricted(φ::P,λ) where P
    _d = φ.d
    _x = φ.z
    λmax = λ
    #=
    x -  λ*d = 0
    λ = x/d
    λmax = minimum(xi/di for i in eachindex(x))
    =#
    for i in 1:50
        break_loop = true
        for i in 1:length(_x)
            if (_x[i] - λmax*_d[i]) < 0
                λmax = 0.5*λmax
                break_loop = false
                break
            end
            break_loop && break
        end
    end
    _x = φ.z
    return λmax
end

"""
    x_sol(res::NLSolvers.ConvergenceInfo)
    
Returns the scalar or vector x that solves the system of equations or is the minimizer of an optimization procedure.
"""
x_sol(res) = NLSolvers.solution(res)
function x_sol(res::NLSolvers.ConvergenceInfo{NLSolvers.BrentMin{Float64}})
    return res.info.x
end

function PositiveLS()
    RestrictedLineSearch(ls_restricted,NLSolvers.Static(0.5))
end
    #ls_restricted(x1,x2) = x2
