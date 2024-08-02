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
    n = loading(model::IsothermModel, p, T)

Calculate the loading based on the model, pressure (p), and temperature (T).

## Inputs
 - model::IsothermModel: the isotherm model
 - p: pressure.
 - T: temperature.

"""
function loading(model::IsothermModel, p, T)
    return loading_ad(model,p,T)
end

function loading_ad(model,p,T)
    return p*ForwardDiff.derivative(p -> p*sp_res(model, p, T), p)
end

"""
    Π = sp_res(model::IsothermModel, p, T)

Calculate the reduced spreading pressure based on the model, pressure (p), and temperature (T). This is also known as the reduced grand potential, and it is defined as:

```math
H = \\int_{0}^{p}\\frac{loading(model,p,T)}{p}dP
```

## Inputs
 - model::IsothermModel: the isotherm model
 - p: pressure.
 - T: temperature.
"""
function sp_res(model, p, T)
    return sp_res_numerical(model, p, T)
end

function sp_res_numerical(model, p, T; solver = QuadGKJL(), abstol = 1e-6, reltol = 1e-6)
        #For cases where the sp_res is not analytical, we use numerical integration

        #Part 1 integral
        ϵ = sqrt(eps(Base.promote_eltype(model,p,T)))

        ∫ni_p⁻¹ = henry_coefficient(model, T)*ϵ

        #Part 2 integral
        f(p) = loading(model, p, T)/p

        prob = IntegralProblem(f(p), (ϵ, p))

        π_i = ∫ni_p⁻¹ + Integrals.solve(prob, solver; reltol = reltol, abstol = abstol)

    return π_i
end

#henry coefficient

#=
loading(model,p) ≈ k*p
dloading/dp* p(p - p0)

=#

"""
    henry_coefficient(model, T)

default units: `[mol/kg]`

Returns the single component spreading pressure of the `model` given the temperature `T`.

## Inputs
 - model::IsothermModel: the isotherm model
 - T: temperature.
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

function pressure_x0(model::IsothermModel, x, T,::typeof(sp_res))
    Π/henry_coefficient(model, T)
end

function pressure_impl(model::IsothermModel, Π, T, ::typeof(sp_res), approx)
    if approx == :exact
        p0 = pressure_x0(model, Π, T)
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
    isosteric_heat(model::IsothermModel, Vᵍ, Vᵃ = zero(eltype(model)), p, T)

Calculate the isosteric heat of adsorption for a given isotherm model.

# Arguments
- `model::IsothermModel`: The isotherm model used to describe the adsorption process.
- `Vᵍ`: The molar volume of the gas phase.
- `Vᵃ`: The molar volume of the adsorbed phase (normaly Vᵃ << Vᵍ, default is zero).
- `p`: Pressure at which the isosteric heat is evaluated.
- `T`: Temperature at which the isosteric heat is evaluated.

# Returns
- The estimated isosteric heat of adsorption.

# Details
The function Estimates the isosteric heat of adsorption for a single component from it's isotherm 
using the Clausius-Clapeyron Equation:

Q_st = T × (Vᵍ - Vᵃ) × (∂n∂T)ₚ/(∂n∂P)ₜ (for explicit loading expressions)

Pan et al. (1998) https://doi.org/10.1021/la9803373

"""
function isosteric_heat(model::IsothermModel, Vᵍ, p, T; Vᵃ = zero(eltype(model)))

    f(∂p,∂T) = loading(model, ∂p, ∂T)
    
    _f,_df = fgradf2(f, p, T)

    ∂n_∂p, ∂n_∂T = _df

    return -T*(Vᵍ - Vᵃ)*∂n_∂T/∂n_∂p
end
