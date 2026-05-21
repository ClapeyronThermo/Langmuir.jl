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
    x = get_adsorbed_composition(iast_solution)
    q_tot = get_q_tot(iast_solution)
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

    #using the last point as initial point for the pressure solver
    #drastically improves speed in models that require integration for sp_res
    Π_old = Π - (sum(x) - 1)*q_tot
    for i in 1:length(isotherms)
        model = isotherms[i]
        if iters == 1
            Pᵢ[i] = pressure(model, Π, T, sp_res)
        else
            p0 = Pᵢ[i]
            Pᵢ[i] = pressure(model, Π, T, sp_res; x0 = Π_old, p0 = p0)
        end
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
