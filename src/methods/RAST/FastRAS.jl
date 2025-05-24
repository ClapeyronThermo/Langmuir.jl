struct FastRAS{IAST<:IASTSolver} <: RASTSolver
    iast::IAST
end

FastRAS(;iast = IASTNestedLoop()) = FastRAS(iast)

function CommonSolve.init(prob::ASTProblem{M,P,TT,Y,G},alg::FastRAS;maxiters = 100,reltol = 1e-12, abstol = 1e-10) where {M,P,TT,Y,G}
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
    conditions = (maxiters,reltol,abstol)
    γ = activity_coefficient(prob.models, prob.T, x)
    iters = 0
    η = similar(x)
    Res = similar(x)
    δ = similar(x)
    Diag = similar(x)
    K = 1 ./ Pᵢ
    η .= 1.0
    state = (;Π,η,K,Diag,Res,δ,γ,x,q_tot,iters,converged)
    return IASTIteration(alg,prob,state,conditions)
end

#=
FastRAS (tentative name? it does not appear on scholar) 
i don't know if this is done before.
in theory, γᵢ depends on ηᵢ, but, we can just squeeze that activity coefficient model in the K value:

RAST equations in fast ias format:

Resᵢ = Πi(ηi) - Πi(η_nc) for i in 1:(nc - 1)
p0ᵢ = ηᵢ/Kᵢ
Res_nc = 1 - sum(Kᵢ*yᵢ*p/ηᵢγᵢ) = 1 - sum(yᵢp/p0ᵢγᵢ)
Diag

=#
function ast_step!(::FastRAS, model::MultiComponentIsothermModel, p, T, y, state::S, maxiters, reltol, abstol) where S
    (;Π,η,K,Diag,Res,δ,γ,x,q_tot,iters,converged) = state
    iters += 1
    #Kpi = scaling factor, p0i = η[i]/K[i]
    n = length(η)
    models = model.isotherms
    ΔJac_nc_nc = zero(eltype(η))
    ΔRes_nc = zero(eltype(η))
    ∑KpiPyiηi = zero(eltype(η))
    Jac_row_last = zero(eltype(η))
    Π_nc = sp_res(last(models), η[end]/K[end], T)
    q⁻¹ = zero(q_tot)
    #update xi
    for i in 1:n
        model = models[i]
        ηᵢ,Kpiᵢ,yᵢ = η[i],K[i],y[i]
        p0ᵢ = ηᵢ/Kpiᵢ
        x[i] = p*yᵢ/p0ᵢ/γ[i]
    end
    γ = activity_coefficient(model, T, x)
    for i in 1:n
        model = models[i]
        ηᵢ,Kpiᵢ,yᵢ = η[i],K[i],y[i]
        p0ᵢ = ηᵢ/Kpiᵢ
        ηᵢ2 = ηᵢ*ηᵢ
        #update last row
        KpiᵢPyᵢ = Kpiᵢ*p*yᵢ/γ[i]
        Jac_rowᵢ = KpiᵢPyᵢ/ηᵢ2
        ∑KpiPyiηi += KpiᵢPyᵢ/ηᵢ
        #Jac_row[i] = Jac_rowᵢ
        qi = loading(model,p0ᵢ,T)
        q⁻¹ += p*yᵢ/qi/p0ᵢ
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
    q_tot = 1/q⁻¹
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
    return (;Π,η,K,Diag,Res,δ,γ,x,q_tot,iters,converged)
end

