"""
    `Sips(M, K₀, E, f)`

    Sips <: IsothermModel

Sips(M, K₀, E, f) represents the Sips isotherm model, which describes the adsorption of a gas on a solid surface.

## Inputs

- `M`::T: maximum loading capacity of the adsorbent, `[mol/kg]`
- `K₀`::T: equilibrium constant at zero coverage, `[1/Pa]`
- `E`::T: adsorption energy, `[J/mol]`
- `f`::T: parameter characterising the heterogeneity of the system  (no units)

## Description

The Sips equation is given by:

n = M * (K₀ * p)^f / (1 + (K₀ * p)^f)

where:
- n is the loading of the adsorbate on the adsorbent,
- M is the maximum loading capacity of the adsorbent,
- K₀ is the equilibrium constant at zero coverage,
- p is the pressure of the gas.

The adsorption energy E is related to the equilibrium constant K₀ by the equation:

K₀ = exp(-E / (R * T))

where:
- R is the gas constant,
- T is the temperature.

"""
struct Sips{T} <: IsothermModel{T}
    M::T
    K₀::T
    E::T
    f::T
end


function sp_res(model::Sips, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    f = model.f
    K = K₀*exp(-E/(Rgas(model)*T))
    return M*log1p((K*p)^f)/f
end

function loading(model::Sips, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    Kpf = (K*p)^f
    return M*Kpf/(1 + Kpf)
end

saturated_loading(model::Sips, T) = model.M #Some depend on T, some don't

function pressure_impl(model::Sips, Π, T,::typeof(sp_res), approx) 
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
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
    l0,p0 = data.l,data.p
    #remove nonzero values
    idx = findall(>(0),l)
    l,p = l0[idx],p0[idx]
    M = maximum(l)
    logp = log.(p)
    loglml = log(l ./ (M .- l))
    _1 = one.(p)
    flogk,f = hcat(_1,logp)\loglml
    logk = flogk/f
    K = exp(logk)
    return T(M,K,zero(K),f)
end

export Sips
