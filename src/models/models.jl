abstract type IsothermModel{T} end

Base.eltype(model::IsothermModel{T}) where T = T
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
model_length(model::IsothermModel) = _model_length(typeof(model))

function _model_length(model::Type{T}) where T <: IsothermModel
    return length(fieldcount(T))
end

"""
    loading(model::IsothermModel, p, T)

Calculate the loading based on the model, pressure (p), and temperature (T).

## Inputs
 - model::IsothermModel: the isotherm model
"""
function loading(model::IsothermModel, p, T)
    return p*ForwardDiff.derivative(p -> p*sp_res(model, p, T), p)
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
    return ForwardDiff.derivative(p -> loading(model, p , T), sqrt(eps(eltype(model))))
end

#inverse problem

"""
sp_res_pressure(model::IsothermModel,q)

given an isotherm::IsothermModel and Π = sp_res(model,p), find p such that sp_res(model,p) = Π.
by default, it performs a root-finding over the isotherm

"""
function sp_res_pressure(model::IsothermModel, Π, T)

    return sp_res_pressure_impl(model, Π, T)
end


function sp_res_pressure_x0(model::IsothermModel, Π, T)
    Π/henry_coefficient(model, T) 
end

function sp_res_pressure_impl(model::IsothermModel, Π, T)
    p0 = sp_res_pressure_x0(model, Π, T)
    f0(p) = Π - sp_res(model, p, T)  
    prob = Roots.ZeroProblem(f0, p0)
    return Roots.solve(prob)
end

function sp_res_pressure_fdf(model, q, T)
    p = sp_res_pressure(model, q, T)
end

#high pressure loading (Is it necessary?)

function saturated_loading(model::IsothermModel, T)
    return loading(model, one(eltype(model))/sqrt(eltype(model)), T)
end

function from_vec(::Type{T},p::AbstractVector{K}) where {T <: IsothermModel,K}
    return T(ntuple(i -> p[i], model_length(model)))
end

function to_vec!(model::IsothermModel,p)
    for i in 1:model_length(model)
        p[i] = getfield(model,i)
    end
    return p
end

export loading, henry_coefficient

include("langmuir.jl")
include("quadratic.jl")
include("bet.jl")
include("henry.jl")
include("temkin.jl")
include("interpolation.jl")

