"""
    `Sips(M, K₀, E, f₀, β)`

    Sips <: IsothermModel

## Inputs

- `M::T`: Saturation loading, `[mol/kg]`
- `K₀::T`: Affinity parameter at T → ∞, `[Pa⁻¹]`
- `E::T`: Adsorption energy, `[J/mol]`
- `f₀::T`: Surface heterogeneity parameter at T → ∞, `[-]`
- `β::T`: Surface heterogeneity coefficient, `[K]`

## Description

The Sips equation is given by:

n = M * (K * p)ᶠ / (1 + (K * p)ᶠ)

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K = K₀ × exp(-E / (R * T))

The exponent f is also temperature dependent and can be expressed as: 

f = f₀ - β/T

where:
- R is the gas constant,
- T is the temperature.

"""
@with_metadata struct Sips{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "affinity parameter")
    (E::T, (-Inf, 0.0), "energy parameter")
    (f₀::T, (0.0, Inf), "surface heterogeneity parameter at T → ∞")
    (β::T, (-Inf, Inf), "surface heterogeneity coefficient")
end


function sp_res(model::Sips, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    f = model.f₀ - model.β/T
    K = K₀*exp(-E/(Rgas(model)*T))
    return M*log1p((K*p)^f)/f
end

function loading(model::Sips, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    f = model.f₀ - model.β/T
    Kpf = (K*p)^f
    return M*Kpf/(1 + Kpf)
end

saturated_loading(model::Sips, T) = model.M #Some depend on T, some don't

function pressure_impl(model::Sips, Π, T,::typeof(sp_res)) 
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    f = model.f₀ - model.β/T
    return expm1(Π*f/M)^(1/f)/K
end

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: Sips
   #l = M*(kp)^f/(1 +(kp)^f)
    #l*(1 +(kp)^f) = M*(kp)^f
    #l + l*(kp)^f = M*(kp)^f
    #M*(kp)^f - l*(kp)^f = l
    #(k*p)^f  = l/(M - l)
    #f*log(k*p) = log(l) - log(M - l)
    #f*(logk) + f*log(p) = log(l/(M - l))
    
    #TODO: f*log(k) = -f*E/RT * log(k0), try to solve for E and K0.
    #remove nonzero values

    Ts, l_p = split_data_by_temperature(data)
    Ms = Vector{eltype(Ts)}(undef, length(l_p))
    logKs = Vector{eltype(Ts)}(undef, length(l_p))
    fs = Vector{eltype(Ts)}(undef, length(l_p))

    for i in eachindex(l_p)
        l_i, p_i = l_p[i]
        idx = findall(>(0.0), l_i)
        l_i, p_i = l_i[idx], p_i[idx]
        M = maximum(l_i) + eps(maximum(l_i)*1.1)
        logp = log.(p_i)
        loglml = log.(l_i ./ (M .- l_i))
        _1 = one.(p_i)
        flogk,f = hcat(_1, logp)\loglml
        logk = flogk/f
        logKs[i] = logk
        Ms[i] = M
        fs[i] = f
    end

    _1 = one(eltype(Ts))
    _1s = ones(eltype(Ts), length(Ts))

    if length(l_p) > 1
        logK0, E = hcat(_1s, _1./ (Rgas(T).*Ts)) \ logKs
        f₀, β = hcat(_1s, -_1./Ts) \ fs
        M = sum(Ms)/length(Ms)
        K0 = exp(logK0)
    else
        M = first(Ms)
        K0 = exp(first(logKs))
        f₀ = first(fs)
        β = 0.0
    end


    return T(M, K0, -E, f₀, β)
end

export Sips
