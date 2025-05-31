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

function loading_at_T(model::IsothermModel, p, T)

    return map(p-> loading(model, p, T),  p)
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

function sp_res_numerical(model::IsothermModel, p, T; solver = QuadGKJL(), abstol = 1e-3, reltol = 1e-3)
    #For cases where the sp_res is not analytical, we use numerical integration

    #Part 1 integral
    ϵ = sqrt(eps(Base.promote_eltype(model, p, T)))

    ∫₁ni_p⁻¹dp = henry_coefficient(model, T)*ϵ

    #Part 2 integral    
    f = let model = model
        (p, T) ->loading(model, p, T)/p
    end

    prob = IntegralProblem(IntegralFunction(f), (ϵ, p), T)

    ∫₂ni_p⁻¹dp = Integrals.solve(prob, solver; reltol = reltol, abstol = abstol).u

    π_i = ∫₁ni_p⁻¹dp + ∫₂ni_p⁻¹dp

return π_i

end

function ChainRulesCore.frule(
    (_, Δmodel, Δp, ΔT), 
    ::typeof(sp_res_numerical), 
    model, 
    p, 
    T; 
    kwargs...
)
return sp_res_numerical(model, p, T, kwargs...), loading(model, p, T)/p
end

@ForwardDiff_frule sp_res_numerical(model, p::ForwardDiff.Dual, T; kwargs...)

#= function sp_res_numerical(model::IsothermModel, p::ForwardDiff.Dual, T; kwargs...)
    return loading(model, p, T)/p
end
 =#

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

pressure_impl(model, x, T) = pressure_impl(model, x, T, sp_res, :exact) # Default to sp_res if f is not provided

# Initial guess for pressure given spreading pressure target
function pressure_x0(model::IsothermModel, Π_target, T,::typeof(sp_res))
    H = henry_coefficient(model, T)
    # Π ≈ H * p  => p ≈ Π / H
    val = Π_target/H
    return ifelse(isfinite(val) & (val > zero(val)), val, one(Π_target)) # Fallback to 1.0 if H is problematic
end

# Initial guess for pressure given loading target
function pressure_x0(model::IsothermModel, q_target, T, ::typeof(loading))
    H = henry_coefficient(model, T)
    # q ≈ H * p => p ≈ q / H
    val = q_target / H
    # If q_target is very small, p_guess might be very small. Ensure it's positive.
    # If H is zero or Inf, this guess can be problematic.
    return ifelse(isfinite(val) & (val > zero(val)), val, one(q_target)) # Fallback to 1.0
end

function pressure_impl(model::IsothermModel, target_val, T, ::typeof(sp_res), approx)
    if approx == :exact
        p0 = pressure_x0(model, target_val, T, sp_res)
        f0 = let model = model, target_val = target_val, T =T
            p -> target_val - sp_res(model, p, T)
        end
        prob = Roots.ZeroProblem(f0, p0)
        sol_p = Roots.solve(prob, Roots.Secant()) # Consider adding bracketing or more robust solver
        return max(0.0, sol_p) # Ensure pressure is not negative
    elseif approx == :henry
        return target_val/henry_coefficient(model, T)
   # elseif approx == :saturated ?
    else
        _0 = Base.promote_eltype(model, target_val, T)
        return _0/_0 # NaN for unimplemented approx
    end
end

function pressure_impl(model::IsothermModel, q_target, T, ::typeof(loading), approx)
    if approx == :exact
        p0 = pressure_x0(model, q_target, T, loading)

        f_to_solve = let model=model, q_target=q_target, T=T
            p -> loading(model, p, T) - q_target
        end

        # Check if target loading is achievable
        sat_load = saturated_loading(model, T)
        # Allow small tolerance for q_target slightly above sat_load due to numerical precision
        if q_target > sat_load * (1 + sqrt(eps(sat_load)))
            @warn "Target loading $q_target is greater than saturated loading $sat_load at temperature $T. Returning Inf."
            return convert(eltype(p0), Inf)
        elseif q_target < 0
             @warn "Target loading $q_target is negative. Returning 0.0."
             return zero(eltype(p0))
        end

        # Safeguard p0 if it's problematic (e.g. Inf/NaN from Henry coeff)
        if !isfinite(p0) || p0 <= 0
            p0 = one(q_target) # Reset to a generic positive guess
        end

        prob = Roots.ZeroProblem(f_to_solve, p0)
        # Consider trying a few initial guesses or a bracketing solver if Secant fails.
        # Roots.Brent() or Roots.AlefeldPotraShi() with a bracket like [0, some_very_high_p]
        # For now, rely on Secant and ensure result is non-negative.
        sol_p = Roots.solve(prob, Roots.Secant()) # Add error handling or maxiters?
        return max(0.0, sol_p) # Ensure pressure is not negative

    elseif approx == :henry
        H = henry_coefficient(model,T)
        return q_target/H
    else
        _0 = Base.promote_eltype(model, q_target, T)
        return _0/_0 # NaN for unimplemented approx
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
function isosteric_heat(model::IsothermModel, p, T; Vᵃ = zero(eltype(p)), Vᵍ = Rgas(model)*T/p)
    
    f =  let model = model
        (∂p,∂T) -> loading(model, ∂p, ∂T)
    end
    
    _f,_df = fgradf2(f, p, T)

    ∂n_∂p, ∂n_∂T = _df

    return T*(Vᵍ - Vᵃ)*∂n_∂T/∂n_∂p
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

export loading, loading_at_T, sp_res, isosteric_heat
export henry_coefficient, saturated_loading
