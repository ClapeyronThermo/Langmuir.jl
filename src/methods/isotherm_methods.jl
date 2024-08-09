#=
api:

necessary function: 
    sp_res(model::IsothermModel, p, T)
    OR
    loading(model::IsothermModel,p,T)


derived:
- loading(model::IsothermModel,p,T) <-> pressure(model,q,T,loading)
- sp_res(model::IsothermModel, p, T) <-> pressure(model,q,T,sp_res)
- henry_coefficient(model::IsothermModel,T), returns Inf if fails.
- saturated_loading(model::IsothermModel,T), returns Inf by default.
=#

"""
    loading(model::IsothermModel, p, T) -> q

Calculate the loading `q` based on the provided isotherm model, pressure `p`, and temperature `T`.

# Arguments
- `model::IsothermModel`: An instance of `IsothermModel`, representing the isotherm model to be used for the calculation.
- `p`: The pressure at which the loading is to be calculated.
- `T`: The temperature at which the loading is to be calculated.

# Returns
- `q`: The calculated loading based on the isotherm model, pressure, and temperature.

# Description
This function computes the loading `q` based on the given isotherm model, pressure `p`,
and temperature `T`.

"""
function loading(model::IsothermModel, p, T)
    return loading_ad(model,p,T)
end

function loading_ad(model,p,T)
    return p*ForwardDiff.derivative(p -> sp_res(model, p, T), p)
end

"""
    sp_res(model::IsothermModel, p, T) -> Π

Calculate the reduced spreading pressure for a given isotherm model at a specific pressure `p` and temperature `T`.

# Arguments
- `model::IsothermModel`: An instance of `IsothermModel`, representing the isotherm model used for the calculation.
- `p`: The pressure at which the reduced spreading pressure is to be calculated.
- `T`: The temperature at which the reduced spreading pressure is to be calculated.

# Returns
- `Π`: The reduced spreading pressure 

# Description
The reduced spreading pressure is a key quantity in Ideal Adsorbed Solution Theory (IAST), used to describe the adsorption behavior of mixtures. This function calculates the reduced spreading pressure Π by integrating the isotherm equation over the pressure range from 0 to `p`.

The reduced spreading pressure is often calculated numerically as:

Π = ∫ (q(p') / p') dp' from 0 to p

where:
- `q(p')` is the loading at pressure `p'`.

"""
function sp_res(model, p, T)
    return sp_res_numerical(model, p, T)
end

function sp_res_numerical(model, p, T; solver = QuadGKJL(), abstol = 1e-6, reltol = 1e-6)
    #For cases where the sp_res is not analytical, we use numerical integration

    #Part 1 integral
    ϵ = sqrt(eps(Base.promote_eltype(model, p, T)))

    ∫₁ni_p⁻¹ = henry_coefficient(model, T)*ϵ

    #Part 2 integral    
    f(p) = loading(model, p, T)/p

    prob = IntegralProblem((u, p) -> f(u), (ϵ, p))

    ∫₂ni_p⁻¹ = Integrals.solve(prob, solver; reltol = reltol, abstol = abstol).u

    π_i = ∫₁ni_p⁻¹ + ∫₂ni_p⁻¹

return π_i

end

#henry coefficient

#=
loading(model,p) ≈ k*p
dloading/dp* p(p - p0)

=#

"""
    henry_coefficient(model::IsothermModel, T) -> H

Calculate the Henry's coefficient for a single component system using the specified isotherm model and temperature `T`.

# Arguments
- `model::IsothermModel`: An instance of `IsothermModel`, representing the isotherm model to be used for the calculation.
- `T`: The temperature at which the Henry's coefficient is to be calculated.

# Returns
- `H`: The Henry's coefficient in the default units of [mol/kg].

# Description
This function returns the Henry's coefficient, which is a measure of the initial slope of the adsorption isotherm at low pressures. It is defined as the derivative of the loading `q` with respect to pressure `p` at `p = 0`:

H = (∂q/∂p) at p = 0 at a given T.
"""
function henry_coefficient(model::IsothermModel, T)
    _0 = zero(eltype(model))
    
    return ForwardDiff.derivative(p -> loading(model, p , T),  _0)
end

#high pressure loading (Is it necessary?)

"""
    saturated_loading(model::IsothermModel, T)

Returns the loading of of an adsorption isotherm when the pressure tends to infinity. by default it is evaluated at `1/√eps(eltype(model))` (`6.7108864e7` for `Float64` inputs.).

## Inputs
 - model::IsothermModel: the isotherm model
 - T: temperature.
"""
function saturated_loading(model::IsothermModel, T)
    ε = eps(Base.promote_eltype(model,T))
    return loading(model, one(ε)/sqrt(ε), T)
end

#inverse problems

"""
pressure(model::IsothermModel, x, T, f; approx = :exact)

given an isotherm::IsothermModel and `x` = f(model,p,T), find `p` such that `f(model,p,T) = x`.
There are two options for `f` and `x`:
- when `f` is `sp_res`, then x = sp_res(model, p, T)
- when `f` is `loading`, then x = loading(model, p, T)

By default, it performs a root-finding over the isotherm. but custom implementations can be done by overloading `pressure_impl(model::MyModel,x,T,f::typeof(f),approx)`
The `approx::Symbol` argument indicates if the procedure is exact or approximate. by default a henry coefficient aproximation is used when `approx =:henry` is used.
"""
function pressure(model::IsothermModel, x, T, f;approx = :exact)
    return pressure_impl(model, x, T, f, approx)
end

pressure_impl(model, x, T) = pressure_impl(model, x, T, sp_res, :exact)

function pressure_x0(model::IsothermModel, Π, T,::typeof(sp_res))
    Π/henry_coefficient(model, T)
end

function pressure_impl(model::IsothermModel, Π, T, ::typeof(sp_res), approx)
    if approx == :exact
        p0 = pressure_x0(model, Π, T, sp_res)
        f0(p,Π) = Π - sp_res(model, p, T)
        prob = Roots.ZeroProblem(f0, p0)
        return Roots.solve(prob,p = Π)
    elseif approx == :henry
        return Π/henry_coefficient(model, T)
   # elseif approx == :saturated ?
    else
        _0 = Base.promote_eltype(model,Π,T)
        return _0/_0
    end
end

"""
    isosteric_heat(model::IsothermModel, Vᵍ, p, T; Vᵃ = zero(eltype(model))) -> Qₛₜ

Calculate the isosteric heat of adsorption for a given isotherm model.

# Arguments
- `model::IsothermModel`: The isotherm model used to describe the adsorption process.
- `Vᵍ`: The molar volume of the gas phase.
- `Vᵃ`: The molar volume of the adsorbed phase (typically Vᵃ << Vᵍ; default is zero).
- `p`: Pressure at which the isosteric heat is evaluated.
- `T`: Temperature at which the isosteric heat is evaluated.

# Returns
- `Qₛₜ`: The estimated isosteric heat of adsorption.

# Description

The function estimates the isosteric heat of adsorption Qₛₜ for a single component using its isotherm and the Clausius-Clapeyron equation:

Qₛₜ = -T * (Vᵍ - Vᵃ) * (∂n/∂T)ₚ / (∂n/∂p)ₜ

where:
- n is the loading,
- Vᵍ is the molar volume of the gas phase,
- Vᵃ is the molar volume of the adsorbed phase,
- T is the temperature,
- p is the pressure.

This equation is derived based on the Clausius-Clapeyron relation, which relates the temperature dependence of the loading to the isosteric heat.

## References:
1. Pan, H., Ritter, J. A., & Balbuena, P. B. (1998). Examination of the approximations used in determining the isosteric heat of adsorption from the Clausius−Clapeyron equation. Langmuir: The ACS Journal of Surfaces and Colloids, 14(21), 6323–6327. [doi:10.1021/la9803373](https://doi.org/10.1021/la9803373)
"""
function isosteric_heat(model::IsothermModel, Vᵍ, p, T; Vᵃ = zero(eltype(model)))
    f(∂p,∂T) = loading(model, ∂p, ∂T)
    
    _f,_df = fgradf2(f, p, T)

    ∂n_∂p, ∂n_∂T = _df

    return -T*(Vᵍ - Vᵃ)*∂n_∂T/∂n_∂p
end

#useful for creating pseudo langmuir models for multicomponent adsoption.
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

export loading, sp_res, isosteric_heat
export henry_coefficient, saturated_loading
