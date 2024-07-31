Base.eltype(model::IsothermModel{T}) where T = T
function Base.eltype(::Type{M}) where M <: IsothermModel{T} where T
    return T
end
#=
api:

necessary function: sp_res(model::IsothermModel, p, T)


derived:
- sp_res_inv(model::IsothermModel,q)
- loading(model::IsothermModel,p) <-> isotherm_pure_pressure(model,q)
- henry_coefficient(model::IsothermModel,p)
- saturated_loading(model::IsothermModel) #TODO: decide if just return max loading or add a flag if the isotherm does not have max loading
=#

Rgas(model) = 8.31446261815324

#default.
model_length(::Type{T}) where T <: IsothermModel = _model_length(T)
model_length(model::IsothermModel) = model_length(typeof(model))

function _model_length(model::Type{T}) where T <: IsothermModel
    return fieldcount(T)
end

"""
    loading(model::IsothermModel, p, T)

Calculate the loading based on the model, pressure (p), and temperature (T).

## Inputs
 - model::IsothermModel: the isotherm model
"""
function loading(model::IsothermModel, p, T)
    return loading_ad(model,p,T)
end

function loading_ad(model,p,T)
    return p*ForwardDiff.derivative(p -> p*sp_res(model, p, T), p)
end

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

Returns the single component spreading pressure of the `model` given the temperature `T`, defined as:

```math
H = `\\lim_{p \\to 0} \\frac{loading(model, p, T)}{p}
```

"""
function henry_coefficient(model::IsothermModel, T)
    _0 = zero(eltype(model))
    
    return ForwardDiff.derivative(p -> loading(model, p , T),  _0)
end

#inverse problem

"""
sp_res_pressure(model::IsothermModel,Π, T, approx = :exact)

given an isotherm::IsothermModel and Π = sp_res(model,p,T), find p such that sp_res(model,p) = Π.
by default, it performs a root-finding over the isotherm.

"""
function sp_res_pressure(model::IsothermModel, Π, T;approx = :exact)
    return sp_res_pressure_impl(model, Π, T, approx)
end

sp_res_pressure_impl(model, Π, T) = sp_res_pressure_impl(model, Π, T, :exact)

function sp_res_pressure_x0(model::IsothermModel, Π, T)
    Π/henry_coefficient(model, T)
end

function sp_res_pressure_impl(model::IsothermModel, Π, T, approx)
    if approx == :exact
        p0 = sp_res_pressure_x0(model, Π, T)
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

function sp_res_pressure_fdf(model, q, T)
    p = sp_res_pressure(model, q, T)
end

#high pressure loading (Is it necessary?)

function saturated_loading(model::IsothermModel, T)
    return loading(model, one(eltype(model))/sqrt(eltype(model)), T)
end

function from_vec(::Type{M},p::AbstractVector{K}) where {M <: IsothermModel,K}
    return M(ntuple(i -> p[i], model_length(M))...)
end

function from_vec(::Type{M},p::NTuple{N,K}) where {M <: IsothermModel,N,K}
    return M(ntuple(i -> p[i], model_length(M))...)
end

function to_vec!(model::IsothermModel,x)
    for i in 1:model_length(model)
        x[i] = getfield(model,i)
    end
    return x
end

function to_vec(model::IsothermModel)
    x = Vector{eltype(model)}(undef, model_length(model))
    to_vec!(model,x)
    return x
end

function to_tuple(model::IsothermModel)
    return ntuple(i -> getfield(model,i), model_length(model))
end

Base.zero(model::M) where M <: IsothermModel = Base.zero(M)

function Base.zero(model::Type{M}) where M <: IsothermModel{T} where T
    from_vec(M,ntuple(Returns(Base.zero(T)),model_length(model)))
end

#when T is not defined (zero(Langmuir))
function Base.zero(model::Type{M}) where M <: IsothermModel
    from_vec(M,ntuple(Returns(0.0),model_length(model)))
end

function Base.iszero(model::IsothermModel)
    result = true
    for i in 1:model_length(model)
        result &= iszero(getfield(model,i))
    end
    return result
end

function x0_guess_fit(::Type{T}, data) where T <: IsothermModel
    eltype = Base.promote_eltype(T, data)
end

export loading, henry_coefficient

include("langmuir.jl")
include("quadratic.jl")
include("bet.jl")
include("henry.jl")
include("temkin.jl")
include("unilan.jl")
include("Toth.jl")
include("interpolation.jl")

