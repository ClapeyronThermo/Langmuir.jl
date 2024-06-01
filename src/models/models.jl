abstract type IsothermModel{T} end

Base.eltype(model::IsothermModel{T}) where T = T
#=
api: 

necessary function: sp_res(model::IsothermModel,p)


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

function loading(model::IsothermModel,p)
    return p*ForwardDiff.derivative(Fix1(sp_res,model),p)
end

#henry coefficient

#=
loading(model,p) ≈ k*p
dloading/dp* p(p - p0)

=#
function henry_coefficient(model::IsothermModel)
    return ForwardDiff.derivative(Fix1(loading,model),sqrt(eps(eltype(model))))
end

#inverse problem

"""
sp_res_pressure(model::IsothermModel,q)

given an isotherm::IsothermModel and q = sp_res(model,p), find p such that sp_res(model,p) = q.
by default, it performs a root-finding
"""
function sp_res_pressure(model::IsothermModel,q)
    return sp_res_pressure_impl(model,p)
end


function sp_res_pressure_x0(model::IsothermModel,q)
    q/henry_coefficient(model)
end

function sp_res_pressure_impl(model::IsothermModel,q)
    p0 = sp_res_pressure_x0(model)
    f0(p) = loading(model,p) - q
    prob = Roots.ZeroProblem(f0,p0)
    return Roots.solve(prob)
end

#high pressure loading
function saturated_loading(model::IsothermModel,p)
    return loading(model,1/sqrt(eltype(model)))
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

include("langmuir.jl")
include("quadratic.jl")
include("bet")
include("henry.jl")
include("temkin.jl")
include("interpolation.jl")

struct Isotherm