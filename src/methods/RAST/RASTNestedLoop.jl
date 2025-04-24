struct RASTNestedLoop{IAST<:IASTSolver} <: RASTSolver
    iast::IAST
end

RASTNestedLoop(;iast = IASTNestedLoop()) = RASTNestedLoop(iast)

function CommonSolve.init(prob::ASTProblem{M,P,TT,Y,G},alg::RASTNestedLoop;maxiters = 100,reltol = 1e-12, abstol = 1e-10) where {M,P,TT,Y,G}
    #solve to a loose level of convergence
    iast_prob = ASTProblem(prob.models.isotherms, prob.p, prob.T, prob.y, prob.x0, prob.gas_model)
    iast_iteration = CommonSolve.init(iast_prob,alg.iast;maxiters,reltol = sqrt(reltol),abstol = sqrt(abstol))
    iast_solution,status = ast_solve!(iast_iteration)
    status == :success || error("IAST guess convergence failed - current number of iterations is $maxiters, consider increasing to meet tolerances.")
    iast_x0 = iast_solution.iter
    Π = iast_x0.Π
    x = iast_x0.x
    q_tot = iast_x0.q_tot
    Pᵢ = get_P0i(iast_solution)
    converged = false
    iters = 0
    conditions = (maxiters,reltol,abstol)
    state = (;Π,Pᵢ,q_tot,x,iters,converged)
    return IASTIteration(alg,prob,state,conditions)
end

#almost equal to IASTNestedLoop, just adding the activity coefficient in the inner solver.
function ast_step!(::RASTNestedLoop, model::MultiComponentIsothermModel, p, T, y, state::S, maxiters, reltol, abstol) where S
    (;Π,Pᵢ,q_tot,x,iters,converged) = state
    iters += 1
    isotherms = model.isotherms

    for i in 1:length(isotherms)
        Pᵢ[i] = pressure(isotherms[i], Π, T, sp_res)
    end

    γᵢ = activity_coefficient(model, T, x)
    q⁻¹ = zero(q_tot)
    for i in 1:length(isotherms)
        puremodel = isotherms[i]
        Pᵢ⁰ = Pᵢ[i]
        xᵢ = y[i]*p/(Pᵢ⁰*γᵢ[i])
        x[i] = xᵢ
        q⁻¹ += xᵢ/loading(puremodel,Pᵢ⁰,T)
    end
    ∑x = sum(x)
    Δ∑x = 1 - sum(x)
    ΔΠ = (∑x - 1)/q⁻¹
    x ./= ∑x
    iters == maxiters && (converged = true)
    abs(ΔΠ) < min(Π*reltol,abstol) && (converged = true)
    Π = Π + ΔΠ
    q_tot = 1/q⁻¹
    return (;Π,Pᵢ,q_tot,x,iters,converged)
end
