abstract type ASTSolver end
abstract type IASTSolver <: ASTSolver end

struct ASTProblem{M,P,TT,Y,X,G}
    models::M
    p::P
    T::TT
    y::Y
    x0::X
    gas_model::G
end

const IASTProblem = ASTProblem

ASTProblem(models, p, T, y; x0 = nothing,gas_model = nothing) = ASTProblem(models, p, T, y, x0, gas_model)

struct ASTIteration{ALG<:ASTSolver,PROB<:ASTProblem,ITER,COND}
    alg::ALG
    prob::PROB
    iter::ITER
    cond::COND
end

const IASTIteration = ASTIteration

function CommonSolve.step!(iter::I) where I <: ASTIteration{ALG,PROB,ITER,COND} where {ALG,PROB,ITER,COND}
    maxiters,reltol,abstol = iter.cond
    alg,prob,state = iter.alg,iter.prob,iter.iter
    new_state = ast_step!(alg,prob.models,prob.p,prob.T,prob.y,state,maxiters,reltol,abstol)
    return ASTIteration(alg,prob,new_state,iter.cond)
end

ast_step!(iter::ASTIteration) = CommonSolve.step!(iter)

get_q_tot(iter::ASTIteration) = get_q_tot(iter.alg,iter)
get_q_tot(alg::ASTSolver,iter::ASTIteration) = iter.iter.q_tot

get_adsorbed_composition(iter::ASTIteration) = get_adsorbed_composition(iter.alg,iter)
get_adsorbed_composition(alg::ASTSolver,iter::ASTIteration) = iter.iter.x

get_P0i(iter::ASTIteration) = get_adsorbed_composition(iter.alg,iter)
get_P0i(alg::ASTSolver,iter::ASTIteration) = iter.iter.Pᵢ


function ast_solve!(state::ASTIteration)
    maxiters = state.cond[1]
    converged = state.iter.converged
    state.iter.converged && return state,:success
    for i in 1:maxiters
        state = step!(state)
        state.iter.converged && return state,:success
    end
    return state,:maxiters_exceeded
end


function CommonSolve.solve!(x0::ASTIteration)
    state,status = ast_solve!(x0)
    q_tot = get_q_tot(state)
    x = get_adsorbed_composition(state)
    return q_tot,x,status
end

function ast_step! end

"""
    iast(models, p, T, y; method = FastIAS(), gas_model = nothing, x0 = nothing, maxiters = 100, reltol = 1e-12, abstol = 1e-10)

Solve for the adsorption equilibrium of a multicomponent gas mixture using the Ideal Adsorbed Solution Theory (IAST).

# Arguments
- `models`: A collection of isotherm models for each component in the mixture.
- `p`: Total pressure of the gas mixture.
- `T`: Temperature of the system.
- `y`: Mole fractions of the components in the bulk gas phase.

# Keyword Arguments
- `method`: Solver method to use (default: `FastIAS()`).
- `gas_model`: Optional gas model for the bulk phase (default: `nothing`).
- `x0`: Initial guess for the adsorbed phase composition (default: `nothing`).
- `maxiters`: Maximum number of iterations allowed (default: `100`).
- `reltol`: Relative tolerance for convergence (default: `1e-12`).
- `abstol`: Absolute tolerance for convergence (default: `1e-10`).

# Returns
- `q_tot`: Total adsorbed amount.
- `x`: Mole fractions of the components in the adsorbed phase.
- `status::Symbol`: Status of the solver (`:success` or `:maxiters_exceeded`).

# Description
The Ideal Adsorbed Solution Theory (IAST) is based on the assumption that the adsorbed phase behaves as an ideal solution. It uses the following key equations:

1. **Equilibrium Condition**:

   ``math
   x_i \\cdot p_i(\\Pi) = y_i \\cdot p
   ``

   where ``\\(x_i\\)`` and ``\\(y_i\\)`` are the mole fractions of component ``\\(i\\)`` in the adsorbed and bulk phases, respectively, ``\\(p_i(\\Pi)\\)`` is the partial pressure of component ``\\(i\\)`` in the adsorbed phase, and ``\\(\\Pi\\)`` is the spreading pressure.

2. **Mass Balance**:

   ``math
   \\sum_i x_i = 1
   ``

3. **Spreading Pressure Relation**:

   ``math
   \\Pi = \\int_0^{p_i} \\frac{q_i(p)}{p} dp
   ``

   where ``\\(q_i(p)\\)`` is the adsorption isotherm for component ``\\(i\\)``.

The `iast` function solves these equations iteratively to determine the adsorbed phase composition ``\\(x\\)`` and the total adsorbed amount ``\\(q_{tot}\\)``.

# Example
```julia
using Langmuir

# Define isotherm models for two components
isotherm_1 = LangmuirS1(1.913, 6.82e-10, -21_976.40)
isotherm_2 = LangmuirS1(1.913, 1.801e-9, -16_925.01)
models = (isotherm_1, isotherm_2)

# Bulk phase conditions
y = [0.5, 0.5]
p = 101325.0
T = 300.0

# Solve using IAST
q_tot, x, status = iast(models, p, T, y)
println("Total adsorbed amount: \$q_tot")
println("Adsorbed phase composition: \$x")
println("Solver status: \$status")
```
"""
function iast end

include("FastIAS.jl")
include("IASTNestedLoop.jl")

function iast(models,p,T,y,method = FastIAS(),gas_model = nothing;x0 = nothing,maxiters = 100,reltol = 1e-12, abstol = 1e-10)
    prob = ASTProblem(models, p, T, y; x0, gas_model)
    return CommonSolve.solve(prob, method; maxiters, reltol, abstol)
end

function iast(models::MultiComponentIsothermModel,p,T,y,method = FastIAS(),gas_model = nothing;x0 = nothing,maxiters = 100,reltol = 1e-12, abstol = 1e-10)
    return iast(models.isotherms,p,T,y,method,gas_model;x0,maxiters,reltol,abstol)
end

#used by rast solvers


export IASTProblem, FastIAS, IASTNestedLoop, iast
