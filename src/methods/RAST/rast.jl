abstract type RASTSolver <: ASTSolver end

"""
    rast(models,p,T,y,method = FastIAS(),gas_model = nothing;x0 = nothing,maxiters = 100,reltol = 1e-12, abstol = 1e-10)

TODO: docs

returns q_tot,x,convergence_symbol (:success, or :maxiters_exceeded)
"""
function rast end

include("RASTNestedLoop.jl")
include("FastRAS.jl")
include("FullRAS.jl")

function rast(models::MultiComponentIsothermModel,p,T,y,method = FullRAS(),gas_model = nothing;x0 = nothing,maxiters = 100,reltol = 1e-12, abstol = 1e-10)
    prob = ASTProblem(models, p, T, y; x0, gas_model)
    return CommonSolve.solve(prob, method; maxiters, reltol, abstol)
end

#fallback
function rast(models,p,T,y,method = RASTNestedLoop(),gas_model = nothing;x0 = nothing,maxiters = 100,reltol = 1e-12, abstol = 1e-10)
    return iast(models,p,T,y,method.iast,gas_model;x0,maxiters,reltol,abstol)
end

export rast, FastRAS, RASTNestedLoop