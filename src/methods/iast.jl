function Δπ!(Δπ, x, p, models)
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

function iast_x0(models,p)
    #langmuir approximation.
    #if we suppose all models are langmuir, then you can obtain an analytical 
    #solution to the iast problem.
    n = length(p)
    M̃ = sum(saturated_loading.(models)) / n # high pressure loading, mean
    K̃ = henry_coefficient.(models) / M̃  # Henry 
    x_guess = K̃ .* p
    x_guess ./= (1 + dot(K̃, p))
    x_guess ./= sum(x_guess)
    return x_guess[1:end-1]
end

function iast(models,p)
    n = length(p)
    x0 = iast_x0(models,p)
    _Δπ!(Δπ, x) = Δπ!(Δπ, FractionVector(x), p, models)
    res = nlsolve(_Δπ!, x0)
    xsol = x_sol(res)
    x = vcat(xsol,1 - sum(xsol))
    p₀ = p ./ x
    nₜ = 1.0 / sum(x[c] / loading(models[c], p₀[c]) for c = 1:n)
    x .*= nₜ
    return x
end