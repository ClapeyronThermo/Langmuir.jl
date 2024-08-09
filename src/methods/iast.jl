function fast_ias_x0(models, p, T, y)
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

function pseudo_langmuir_params(model, p, T, Πmin, Πmax)
    M = saturated_loading(model, T)
    Kh = henry_coefficient(model, T)
    
    if isfinite(Kh) && isfinite(M)
        K = Kh/M
        #K/Kh = (Kh/M)/Kh = 1/M
        #PKave = P*mean(K*M)/Mi
        return M,K
    else
        pmin = pressure(model, Πmin, T, sp_res)
        pmax = pressure(model, Πmax, T, sp_res)
        lmin = loading(model,pmin,T)
        lmax = loading(model,pmax,T)
        lvec = SVector((lmin,lmax))
        pvec = SVector((pmin,pmax))
        _MK,_K = hcat(pvec,-lvec .* pvec)\lvec
        _M = _MK/_K
        return _M,_K
    end
end

function iast_Π0(models, p, T, y, x0 = nothing)
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

function iast(models, p, T, y; x0 = nothing,ss_iters = 3*length(y), fastias_iters = 100)
    n = length(models)
    #TODO: fastIAS
    Π0,p_i = iast_Π0(models, p, T, y, x0)
    Πx,x = iast_nested_loop(models, p, T, y, Π0, p_i, ss_iters)
    return Πx, x
end

function iast_nested_loop(models::M, p, T, y, Π0, p_i = similar(y), iters = 5) where M
    Π = Π0
    x = similar(p_i)
    for k in 1:iters
        for i in 1:length(p_i)
            model = models[i]
            p_i[i] = pressure(model, Π, T, sp_res)
        end
        q_tot_inv = zero(eltype(p_i))
        for i in 1:length(p_i)
            model = models[i]
            p0i = p_i[i]
            x[i] = p*y[i]/p0i
            qi = loading(model,p0i,T)
            q_tot_inv +=p*y[i]/p0i/qi
        end
        Δ∑x = 1 - sum(x)
        Δ2 = -Δ∑x/q_tot_inv
        if abs(Δ∑x) < sqrt(eps(eltype(Δ∑x)))
            break
        end
        Π = Π + Δ2
    end
    return Π,x
end

function update_ias_loop!(models, p, T, y, η, Res, Jac, δ, Kpi)
    Jac_col,Jac_row,Jac_diag = Jac
    n = length(η)

    ΔJac_nc_nc = zero(eltype(η))
    ΔRes_nc = zero(eltype(η))
    ∑KpiPyiηi = zero(eltype(η))
    Π_nc = sp_res(last(models), η[end],T)
    for i in 1:n
        model = models[i]
        ηi = η[i]
        ηi2 = ηi*ηi
        #update last row
        KpiPyi = Kpi[i]*p*y[i]
        Jac_rowi = KpiPyi/ηi2
        ∑KpiPyiηi += KpiPyi/ηi
        Jac_row[i] = Jac_rowi
        qi = loading(model,ηi,T)
        #update diagonals
        Jac_diagi = qi/ηi
        Jac_diag[i] = Jac_diagi
        if i < n
            Resi = sp_res(model, ηi, T) - Π_nc
            Res[i] = Resi
            ΔRes_nc += Resi*Jac_rowi/Jac_diagi
            ΔJac_nc_nc += Jac_rowi/Jac_diagi
        end
    end

    #update last term of the last row of the jac
    Jac_nc_nc = Jac_row[end] + Jac_diag[end]*ΔJac_nc_nc
    Jac_diag[end] = Jac_nc_nc
    Jac_row[end] = Jac_nc_nc

    #update last term of residual
    Res[end] = 1 - ∑KpiPyiηi - ΔRes_nc


    #solve system of equations by backsubstitution
    δ_nc = -Res[end]/Jac_diag[end]
    δ[end] = δ_nc

    for i in 1:(n-1)
        δ[i] = (-Res[i] + Jac_diag[i]*δ_nc)/Jac_diag[i]
    end

    #update η
    for i in 1:n
        ηi = η[i]
        δi = δ[i]
        if ηi + δi < 0
            η[i] = 0.5*ηi
        else
            η[i] = ηi + δi
        end
    end
    
    J = similar(Jac_diag,(3,3))
    J .= 0
    for i in 1:3
        J[i,i] = Jac_diag[i]
        J[i,end] = Jac_row[i]
    end
    return J,Res,δ
    


    #update diag
end

function fast_ias(models::M,p,T,y,x0 = nothing) where M
    if x0 === nothing
        Kpi,η = fast_ias_x0(models, p, T, y)
    else
        Kpi,η = fast_ias_x0(models, p, T, y, x0)
    end
end

[2.847
0.028
0.
2.223
1.228
0.]

