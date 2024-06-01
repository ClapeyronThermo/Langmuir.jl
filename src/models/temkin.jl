struct TemkinApprox{T} <: IsothermModel
    M::T
    K::T
    theta::T
end

function sp_res(model::TemkinApprox,p)
    M,K,θ = model.M,model.K,model.theta
    Kp = K*p
    return M*(log1p(Kp) + θ*(2*Kp + 1)/(2*(Kp+1)*(Kp+1)))
end

export TemkinApprox