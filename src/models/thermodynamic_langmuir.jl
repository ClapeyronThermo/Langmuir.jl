@with_metadata struct ThermodynamicLangmuir{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "affinity parameter") 
    (E::T, (-Inf, 0.0), "energy parameter")
    (Bᵢᵩ::T, (-Inf, 0.0), "adsorbate-adsorbent interaction coefficient")
end


function gibbs_excess_free_energy(model::ThermodynamicLangmuir, T, θ)

    _1 = one(eltype(T))
    Bᵢᵩ = model.Bᵢᵩ
    T⁻¹ = _1/T
    θᵢ, θᵩ = θ 

    τᵢᵩ = Bᵢᵩ*T⁻¹
    Gᵢᵩ = exp(-0.3*τᵢᵩ)

    gᴱ_RT = θᵢ*θᵩ*τᵢᵩ*(Gᵢᵩ - _1)/(θᵢ*Gᵢᵩ + θᵩ)

    return gᴱ_RT
end

function activity_coefficient(model::ThermodynamicLangmuir, T, θ)

         fun(θ) = gibbs_excess_free_energy(model, T, θ)

         return exp.(gradient(fun, θ))
end

function loading_x0(model::ThermodynamicLangmuir, p, T)    

    guess_model = from_vec(LangmuirS1, [model.M, model.K₀, model.E])

    return loading(guess_model, p, T)
end

function loading_impl(model::ThermodynamicLangmuir, p, T)

    _1 = one(eltype(p))
    q0 = loading_x0(model, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))

    f0(q, P_T) = begin 
    p, T = P_T
    θᵢ = q/M
    γᵢ, γᵩ = activity_coefficient(model, T, [θᵢ, _1 - θᵢ])
    q - M * K *p / (γᵢ/γᵩ + K*p)    
    end 
    prob = Roots.ZeroProblem(f0, q0)
    return Roots.solve(prob, p = (p, T))

end


function loading(model::ThermodynamicLangmuir, p, T)
    return loading_impl(model, p, T)
end

export ThermodynamicLangmuir

