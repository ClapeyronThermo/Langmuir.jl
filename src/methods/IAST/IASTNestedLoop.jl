struct IASTNestedLoop <: IASTSolver end

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
