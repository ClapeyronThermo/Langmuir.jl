function fast_ias_loop!(models,η)
    #x::Vector{Float64}, 
    #p::Vector{Float64},
    #aim::Vector{<:AdsIsoTModel})
    n = length(x)
    @assert length(Δπ) == n - 1
    sp_res_n = sp_res(last(models), p[end]/x[end])
    for i = 1:(n-1)
        Δπ[i] = sp_res(model[i], p[i]/x[i]) - sp_res_n
    end
    return Δπ
end

function iast_Π0(models,p,y,x0)
    if x0 === nothing
        K_henry = henry_coefficient.(models)
        K_average = dot(K_henry,y)
        Π0 = zero(eltype(K_average))*p
        #Mangano et al. 2015, initial guess (eq. 16)
        for i in 1:length(y)
            model_i = models[i]
            p0i = P*K_average/K_henry[i]
            Π0 = min(Π0,sp_res(model_i, p0i))
        end
        return Π0
    else
        #xi*loading(mi,Π0) = py

    end
end

function iast(models,p,y;x0 = nothing,ss_iters = 3*length(y),fastias_iters = 100)
    n = length(p)
    #TODO: fastIAS
    Π0 = iast_Π0(models,p,y,x0)
    p_i = similar(y)
    Πx = iast_nested_loop(models,p,y,Π0,p_i,ss_iters)
    x = similar(y)
    for i in 1:n
        x[i] = y[i]p/p_i[i]
    end
    x[i] ./= sum(x)
    return Πx,x
end

function iast_nested_loop(models::T,p,y,Π,p_i = similar(y),iters = 5) where T
    function iast_f0(Π) 
        f = one(Π)
        df = zero(Π)
        for i in 1:length(y)
            mi = models[i]
            p0i = sp_res_pressure(mi,Π)
            p_i[i] = p0i
            fi = p*y[i]/p0i
            f -= fi
            df -= fi/loading(mi,p0i)
        end
        return f,f/df
    end
    prob = Roots.ZeroProblem(iast_f0,Π)
    #newton with history
    return Roots.solve(prob,Roots.LithBoonkkampIJzerman(4,1),maxiters = iters)
end