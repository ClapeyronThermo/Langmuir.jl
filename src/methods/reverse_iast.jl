function reverse_iast(models,p,x,y0 = nothing;ss_iters = 100)
    Π0 = reverse_iast_Π0(models,p,x,y0)
    p_i = similar(x)
    Πy = reverse_iast_nested_loop(models,p,y,Π0,p_i,ss_iters)
    y = similar(x)
    for i in 1:n
        x[i] = x[i]*p_i[i]/p
    end
    y[i] ./= sum(y)
    return Πy,y
end

function reverse_iast_Π0(models,p,x,y0 = nothing)
    #=
    lets suppose each model is a henry isotherm
    xi*p_i(Π) = yi*p
    xi*Khi/Π = yi*p
    xi*Khi/p = yi*Π
    Π = sum(xi*Khi)/p
    =#
    Kh = henry_coefficient.(models)
    return dot(Kh,x)/p
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

function reverse_iast_nested_loop(models::T,p,x,Π,p_i = similar(y),iters = 5) where T
    function iast_f0(Π) 
        f = one(Π)
        df = zero(Π)
        for i in 1:length(y)
            mi = models[i]
            p0i = sp_res_pressure(mi,Π)
            p_i[i] = p0i
            fi = x[i]*p0i/p
            f -= fi #TODO: verify if this is correct
            df -= fi*loading(mi,p0i)
        end
        return f,f/df
    end
    prob = Roots.ZeroProblem(iast_f0,Π)
    #newton with history
    return Roots.solve(prob,Roots.LithBoonkkampIJzerman(4,1),maxiters = iters)
end