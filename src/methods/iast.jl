abstract type ASTSolver end
abstract type IASTSolver <: ASTSolver end
struct FastIAS <: IASTSolver end
struct IASTNestedLoop <: IASTSolver end

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

get_P0i(iter::ASTIteration{IASTNestedLoop}) = iter.iter.Pᵢ
get_P0i(iter::ASTIteration{FastIAS}) = iter.iter.η ./ iter.iter.K

const IASTIteration = ASTIteration

function CommonSolve.step!(iter::I) where I <: ASTIteration{ALG,PROB,ITER,COND} where {ALG,PROB,ITER,COND}
    maxiters,reltol,abstol = iter.cond
    alg,prob,state = iter.alg,iter.prob,iter.iter
    new_state = ast_step!(alg,prob.models,prob.p,prob.T,prob.y,state,maxiters,reltol,abstol)
    return ASTIteration(alg,prob,new_state,iter.cond)
end

ast_step!(iter::ASTIteration) = CommonSolve.step!(iter)

function ast_solve!(x::ASTIteration)
    maxiters = x.cond[1]
    converged = x.iter.converged
    x.iter.converged && return x,:success
    for i in 1:maxiters
        x = step!(x)
        x.iter.converged && return x,:success
    end
    return x,:maxiters_exceeded
end

function CommonSolve.solve!(x0::ASTIteration)
    x,status = ast_solve!(x0)
    return x.iter.q_tot,x.iter.x,status
end

function iast_x0(::FastIAS, models, p, T, y, x0 = nothing)
    n = length(models)
    #we calculate Πi, find the extrema, and interpolate a langmuir isotherm.
    #this skips the non existent henry and/or saturated_loading isotherms when they are not available
    Πi = sp_res.(models,p,T)
    Πmin, Πmax = extrema(Πi)
    MK = pseudo_langmuir_params.(models, p, T, Πmin, Πmax)
    #return Kpi
    η = similar(y,eltype(Πi))
    K = similar(y,eltype(Πi))
    K .= last.(MK)
    if !isnothing(x0)
        ∑x0 = sum(x0)
        for i in eachindex(η)
            p0i = ∑x0*p*y[i]/x0[i]
            η[i] = p0i*K[i]
        end
        return η,K
    end

    Kh_avg = zero(eltype(η))
    for i in eachindex(y)
        Mᵢ,Kᵢ = MK[i]
        Kh_avg += y[i]*Mᵢ*Kᵢ
    end
    Kh_avg /= sum(y)
    pKh_avg = Kh_avg
    for i in eachindex(η)
        Mᵢ,Kᵢ = MK[i]
        KHᵢ = Kᵢ * Mᵢ
        s1 =  p*Kᵢ
        s2 = p*Kh_avg/Mᵢ
        η[i] = min(s1,s2)
    end
    return η,K
end

function iast_x0(::IASTNestedLoop,models, p, T, y, x0 = nothing)
    if !isnothing(x0)
        P0i = similar(y,Base.promote_eltype(models[1],p,T,y,x0))
        Π0 = Inf*one(eltype(P0i))
        for i in eachindex(y)
            model_i = models[i]
            ∑x0 = sum(x0)
            p0i = ∑x0*p*y[i]/x0[i]
            P0i[i] = p0i
            Π0 = min(Π0,sp_res(model_i, p0i, T))
        end
        return Π0,P0i
    end
    Πi = sp_res.(models,p,T)
    Πmin, Πmax = extrema(Πi)
    MK = pseudo_langmuir_params.(models, p, T, Πmin, Πmax)

    Kh_avg = zero(eltype(Πmin))
    for i in eachindex(y)
        Mᵢ,Kᵢ = MK[i]
        Kh_avg += y[i]*Mᵢ*Kᵢ
    end
    Π0 = Inf*Kh_avg
    P0i = zeros(eltype(Πmin),length(y))
    for i in eachindex(y)
        Mᵢ,Kᵢ = MK[i]
        KHᵢ = Kᵢ * Mᵢ
        model_i = models[i]
        p0i = p*Kh_avg/KHᵢ
        P0i[i] = 1/KHᵢ
        Π0 = min(Π0,sp_res(model_i, p0i, T))
    end
    P0i .*= Π0
    return Π0,P0i
end

function CommonSolve.init(prob::ASTProblem{M,P,TT,Y,G},alg::IASTNestedLoop;maxiters = 100,reltol = 1e-12, abstol = 1e-10) where {M,P,TT,Y,G}
    Π,Pᵢ = iast_x0(alg,prob.models, prob.p, prob.T, prob.y, prob.x0)
    q_tot = zero(eltype(Pᵢ))/zero(eltype(Pᵢ))
    x = similar(Pᵢ)
    iters = 0
    converged = false
    conditions = (maxiters,reltol,abstol)
    state = (;Π,Pᵢ,q_tot,x,iters,converged)
    return IASTIteration(alg,prob,state,conditions)
end

function CommonSolve.init(prob::ASTProblem{M,P,TT,Y,G},alg::FastIAS;maxiters = 100,reltol = 1e-12, abstol = 1e-10) where {M,P,TT,Y,G}
    η,K = iast_x0(alg,prob.models, prob.p, prob.T, prob.y, prob.x0)
    q_tot = zero(eltype(η))/zero(eltype(η))
    x = similar(η)
    Diag = similar(η)
    Res = similar(η)
    δ = similar(η)
    iters = 0
    converged = false
    conditions = (maxiters,reltol,abstol)
    Π = NaN*zero(eltype(η))
    state = (;Π,η,K,Diag,Res,δ,x,q_tot,iters,converged)
    return IASTIteration(alg,prob,state,conditions)
end

function ast_step!(::IASTNestedLoop, models, p, T, y, state::S, maxiters, reltol, abstol) where S
    (;Π,Pᵢ,q_tot,x,iters,converged) = state
    iters += 1

    for i in 1:length(Pᵢ)
        model = models[i]
        Pᵢ[i] = pressure(model, Π, T, sp_res)
    end

    q⁻¹ = zero(eltype(Pᵢ))

    for i in 1:length(Pᵢ)
        model = models[i]
        p0i = Pᵢ[i]
        x[i] = p*y[i]/p0i
        qi = loading(model,p0i,T)
        q⁻¹ +=p*y[i]/p0i/qi
    end

    ΔΠ = (sum(x) - 1)/q⁻¹
    iters == maxiters && (converged = true)
    abs(ΔΠ) < min(Π*reltol,abstol) && (converged = true)
    Π = Π + ΔΠ
    q_tot = 1/q⁻¹
    return (;Π,Pᵢ,q_tot,x,iters,converged)
end

function ast_step!(::FastIAS, models, p, T, y, state::S, maxiters, reltol, abstol) where S
    (;Π,η,K,Diag,Res,δ,x,q_tot,iters,converged) = state
    iters += 1
    #Kpi = scaling factor, p0i = η[i]/K[i]
    n = length(η)
    ΔJac_nc_nc = zero(eltype(η))
    ΔRes_nc = zero(eltype(η))
    ∑KpiPyiηi = zero(eltype(η))
    Jac_row_last = zero(eltype(η))
    Π_nc = sp_res(last(models), η[end]/K[end], T)
    q_tot_inv = zero(q_tot)
    for i in 1:n
        model = models[i]
        ηᵢ,Kpiᵢ,yᵢ = η[i],K[i],y[i]
        p0ᵢ = ηᵢ/Kpiᵢ
        ηᵢ2 = ηᵢ*ηᵢ
        #update last row
        KpiᵢPyᵢ = Kpiᵢ*p*yᵢ
        Jac_rowᵢ = KpiᵢPyᵢ/ηᵢ2
        ∑KpiPyiηi += KpiᵢPyᵢ/ηᵢ
        #Jac_row[i] = Jac_rowᵢ
        qi = loading(model,p0ᵢ,T)
        q_tot_inv += p*yᵢ/qi/p0ᵢ
        x[i] = p*yᵢ/p0ᵢ
        #update diagonals
        Diagᵢ = qi/ηᵢ
        Diag[i] = Diagᵢ
        if i != n
            Resᵢ = sp_res(model,p0ᵢ, T) - Π_nc
            Res[i] = Resᵢ
            ΔRes_nc += Resᵢ*Jac_rowᵢ/Diagᵢ
            ΔJac_nc_nc += Jac_rowᵢ/Diagᵢ
        else
            Jac_row_last = Jac_rowᵢ
        end
    end
    q_tot = 1/q_tot_inv
    #update last term of the last row of the jac
    Jac_nc_nc = Jac_row_last + Diag[end]*ΔJac_nc_nc
    Jac_row_nc = Diag[end]
    Diag[end] = Jac_nc_nc

    #update last term of residual
    Res[end] = 1 - ∑KpiPyiηi - ΔRes_nc

    #solve system of equations by backsubstitution
    δ_nc = -Res[end]/Diag[end]
    δ[end] = δ_nc
    for i in 1:(n-1)
        δ[i] = -(Res[i] - Jac_row_nc*δ_nc)/Diag[i]
    end

    #update η
    norm_η = -Inf*one(eltype(η))
    for i in 1:n
        ηi = η[i]
        δi = δ[i]
        if ηi + δi < 0
            norm_η
            η[i] = 0.5*ηi
            δ[i] = -0.5*δi

        else
            η[i] = ηi + δi
        end
    end
    ΔRes = norm(δ,Inf)
    ΔRes <= abstol && (converged = true)
    norm(δ,1) <= reltol && (converged = true)
    Π = Π_nc
    return (;Π,η,K,Diag,Res,δ,x,q_tot,iters,converged)
end


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
function iast(models,p,T,y,method = FastIAS(),gas_model = nothing;x0 = nothing,maxiters = 100,reltol = 1e-12, abstol = 1e-10)
    prob = ASTProblem(models, p, T, y; x0, gas_model)
    return CommonSolve.solve(prob, method; maxiters, reltol, abstol)
end

function iast(models::ThermodynamicLangmuir,p,T,y,method = FastIAS(),gas_model = nothing;x0 = nothing,maxiters = 100,reltol = 1e-12, abstol = 1e-10)
    return iast(models.isotherms,p,T,y,method,gas_model;x0,maxiters,reltol,abstol)
end


export IASTProblem, FastIAS, IASTNestedLoop, iast
