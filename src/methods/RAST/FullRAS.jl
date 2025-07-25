struct FullRAS{IAST<:IASTSolver} <: RASTSolver
    iast::IAST
end

FullRAS(;iast = IASTNestedLoop()) = FullRAS(iast)

function CommonSolve.init(prob::ASTProblem{M,P,TT,Y,G},alg::FullRAS;maxiters = 100,reltol = 1e-12, abstol = 1e-10) where {M,P,TT,Y,G}
    #solve to a loose level of convergence
    iast_prob = ASTProblem(prob.models.isotherms, prob.p, prob.T, prob.y, prob.x0, prob.gas_model)
    iast_iteration = CommonSolve.init(iast_prob,alg.iast;maxiters,reltol = sqrt(reltol),abstol = sqrt(abstol))
    iast_solution,status = ast_solve!(iast_iteration)
    status == :success || error("IAST guess convergence failed - current number of iterations is $maxiters, consider increasing to meet tolerances.")
    iast_x0 = iast_solution.iter
    Π = iast_x0.Π
    x = get_adsorbed_composition(iast_solution)
    Pᵢ = get_P0i(iast_solution)
    converged = false
    conditions = (maxiters,reltol,abstol)
    
    iters = 0
    η = similar(x)
    Res = similar(x)
    δ = similar(x)
    Diag = similar(x)
    K = 1 ./ Pᵢ
    η .= 1.0
    n = length(η)
    f!(_F,_z) = full_ras_system(_F,_z,prob.models,prob.p,prob.T,prob.y,K)
    F = similar(x,2n)
    J = similar(x,(2n,2n))
    J2 = similar(J)
    z = similar(F)
    s = similar(F)
    piv = zeros(Int,2n)
    config = ForwardDiff.JacobianConfig(f!,F,z)
    cache = (f!,F,J,z,config)
    
    z[1:n] .= η
    z[n+1:end] .= x
    state = (;K,f!,F,J,J2,piv,z,s,config,iters,converged)
    return IASTIteration(alg,prob,state,conditions)
end

#=
FullRAS

we use the FastIAS framework (optimizing in terms of η = p0ᵢ*Kᵢ), but fully incorporating the activity coefficient.

variables: vcat(ηi,xi)
Resᵢ = Πᵢ(ηi) - Πᵢ(η_nc), i in 1:nc - 1
Resᵢ = xᵢ*ηᵢγᵢ - Kᵢ*yᵢ*p, i in nc:(2nc - 1)
Resᵢ = 1 - sum(Kᵢ*yᵢ*p/ηᵢγᵢ), i = 2nc
=#

function full_ras_system(F,z,system::MultiComponentIsothermModel, p, T, y, K)
    models = system.isotherms
    n = length(K)
    η = @view z[1:n]
    x = @view z[(n+1):end]
    γ = activity_coefficient(system, T, x)
    ∑x = one(eltype(F))
    F1 = @view F[1:n]
    F2 = @view F[(n+1):(2n)]
    for i in 1:n
        model = models[i]
        ηᵢ,Kpiᵢ,yᵢ = η[i],K[i],y[i]
        xᵢ,γᵢ = x[i],γ[i]
        p0ᵢ = ηᵢ/Kpiᵢ
        #evaluate sp_res
        Πᵢ = sp_res(model,p0ᵢ,T)
        F1[i] = Πᵢ
        F2[i] = xᵢ*ηᵢ*γᵢ - Kpiᵢ*yᵢ*p
        ∑x -= yᵢ*p/(p0ᵢ*γᵢ)
    end
    Πₙ = F1[n]
    for i in 1:n-1
        F1[i] -= Πₙ
    end
    F1[end] = ∑x
    return F
end

function ast_step!(::FullRAS, model::MultiComponentIsothermModel, p, T, y, state::S, maxiters, reltol, abstol) where S
    (;K,f!,F,J,J2,piv,z,s,config,iters,converged) = state
    n = length(x)
    iters += 1
    #
    ForwardDiff.jacobian!(J,f!,F,z,config,Val{false}())
    J2 .= J
    lu = unsafe_LU!(J,piv)
    s .= -F
    ldiv!(lu,s)
    for i in 1:length(s)
        si = s[i]
        abs(si) < eps(eltype(s)) && si < 0 && (s[i] = 0)
        new_z = z[i] + si
        if new_z < 0
            z[i] = z[i]/2
        else
            z[i] = new_z
        end
    end
    Fnorm = 0.5*dot(F,F)
    converged = Fnorm < abstol
    ΔRes = norm(s,Inf)
    ΔRes <= abstol && (converged = true)
    norm(δ,1) <= reltol && (converged = true)
    return (;K,f!,F,J,J2,piv,z,s,config,iters,converged)
end

function get_P0i(alg::ASTSolver,iter::FullRAS) 
    K = iter.iter.K
    n = length(K)
    z = iter.iter.z
    η = @view z[1:n]
    return η ./ K
end

function get_adsorbed_composition(alg::ASTSolver,iter::FullRAS) 
    z = iter.iter.z
    return z[n+1:end]
end

function get_q_tot(alg::ASTSolver,iter::FullRAS) 
    prob = iter.prob
    models = prob.models.isotherms
    K = iter.iter.K
    z = iter.iter.z
    q⁻¹ = zero(eltype(z))
    for i in 1:length(models)
        puremodel = models[i]
        Pᵢ⁰ = z[i]/K[i]
        q⁻¹ += x[i]/loading(puremodel,Pᵢ⁰,T)
    end
    return 1/q⁻¹
end

export FullRAS