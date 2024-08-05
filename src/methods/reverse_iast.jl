function reverse_iast(models,p,T,x,Π0 = reverse_iast_Π0(models,p,T,x);ss_iters = 100)
    #Π0 = reverse_iast_Π0(models,p,T,x,y0)
    p_i = similar(x)
    @show Π0
    Πy = reverse_iast_nested_loop(models,p,T,x,Π0,p_i,ss_iters)
    y = similar(x)
    for i in 1:length(models)
        y[i] = y[i]*p_i[i]/p
    end
    y ./= sum(y)
    return Πy,y
end

function reverse_iast_Π0(models,p,T,x,y0 = nothing)
    #=
    lets suppose each model is a henry isotherm
    xi*p_i(Π) = yi*p
    xi*Khi/Π = yi*p
    xi*Khi/p = yi*Π
    Π = sum(xi*Khi)/p
    =#
    Kh = henry_coefficient.(models,T)
    return minimum(Kh)*p
    return dot(Kh,x)*p
end
#nested loop approach:
#=1 - sum(yi) = 0
iast conditions: xi*p_i(Π) = yi*p
yi = xi*p_i(Π)/p
f(Π) = 1 - sum(xi*p_i(Π)/p) = 0
df(Π)/dΠ = -sum(xi*pi(Π))
p_i = p_i - (Π(p_i) - Πset)/dΠ/dp
dpi/dΠ = 1/dΠ/dp


=#

function reverse_iast_nested_loop(models::TT,p,T,x,Π,p_i = similar(y),iters = 5) where TT
    function iast_f0(logΠ)
        Π = exp(logΠ) 
        f = one(Π)
        df = zero(Π)
        for i in 1:length(x)
            mi = models[i]
            p0i = sp_res_pressure(mi,Π,T)
            p_i[i] = p0i
            fi = x[i]*p0i/p
            f -= fi #TODO: verify if this is correct
            df -= fi*loading(mi,p0i,T)
            @show f,df
        end
        return f,f/df
    end
    prob = Roots.ZeroProblem(iast_f0,log(Π))
    logΠ = Roots.solve(prob, Roots.LithBoonkkampIJzerman(4,1), maxiters = iters)
    return exp(logΠ)
end