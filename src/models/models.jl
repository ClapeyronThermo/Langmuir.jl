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
    return loading_ad(model, p, T)
end

function loading_ad(model, p, T)
    return p*ForwardDiff.derivative(p -> sp_res(model, p, T), p)
end

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

"""
    sp_res_pressure(model::IsothermModel, Π, T) -> p

Find the pressure `p` such that `sp_res(model, p)` equals the given residual pressure `Π`.

"""
function sp_res_pressure(model::IsothermModel, Π, T)
    return sp_res_pressure_impl(model, Π, T)
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

function saturated_loading(model::IsothermModel, T)
    return loading(model, one(eltype(model))/sqrt(eltype(model)), T)
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

include("langmuir.jl")
include("quadratic.jl")
include("bet.jl")
include("henry.jl")
include("temkin.jl")
include("unilan.jl")
include("Toth.jl")
include("interpolation.jl")

