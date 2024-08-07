Base.eltype(model::IsothermModel{T}) where T = T

function Base.eltype(::Type{M}) where M <: IsothermModel{T} where T
    return T
end

function (::Type{I})(p) where I <: IsothermModel
    return from_vec(I,p)
end

Rgas(model) = 8.31446261815324

#default.
model_length(::Type{T}) where T <: IsothermModel = _model_length(T)
model_length(model::IsothermModel) = model_length(typeof(model))

function _model_length(model::Type{T}) where T <: IsothermModel
    return fieldcount(T)
end

function from_vec(m::IsothermModel,x)
    return from_vec(typeof(m),x)
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

Reference:
- Pan et al. (1998) DOI: 10.1021/la9803373
"""
function isosteric_heat(model::IsothermModel, Vᵍ, p, T; Vᵃ = zero(eltype(model)))

        f(∂p,∂T) = loading(model, ∂p, ∂T)
        
        _f,_df = fgradf2(f, p, T)

        ∂n_∂p, ∂n_∂T = _df

        return -T*(Vᵍ - Vᵃ)*∂n_∂T/∂n_∂p
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

export loading, henry_coefficient, isosteric_heat, sp_res

include("freundlich.jl")
include("langmuir.jl")
include("langmuir_freundlich.jl")
include("redlich_peterson.jl")
include("sips.jl")
include("quadratic.jl")
include("bet.jl")
include("henry.jl")
include("temkin.jl")
include("unilan.jl")
include("toth.jl")
include("multisite.jl")
include("interpolation.jl")
