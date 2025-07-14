"""
    ThermodynamicLangmuir(M, K₀, E, Bᵢᵩ)

The `ThermodynamicLangmuir` struct represents a thermodynamic Langmuir model with activity coefficients calculated using an aNRTL-like approach.

## Parameters

- `M`: Saturation loading `[mol/kg]`
- `K₀`: Affinity parameter at high temperature `[1/Pa]`
- `E`: Adsorption energy `[J/mol]`
- `Bᵢᵩ`: Interaction energy parameter between adsorbate species `i` and vacant sites `[J/mol]`

## Description

The Langmuir equation is given by:

    nᵢ = (M * K * P) / (γᵢ/γᵩ + K * P)

where `nᵢ` is the adsorbed amount of species `i`, `K` is calculated as:

    K = K₀ * exp(-E / (R * T))

where:
- `R` is the gas constant
- `T` is the temperature

The activity coefficients `γᵢ` and `γᵩ` are determined using the Gibbs excess free energy, `gᴱ/RT`, which is calculated based on the surface fractions (`θᵢ`, `θᵩ`) and interaction parameters derived from `Bᵢᵩ`. This free energy value is used in the `activity_coefficient` function to compute the activity coefficients of the adsorbate and phantom molecules.
"""
@with_metadata struct ThermodynamicLangmuir{T} <: IsothermModel{T}
    (M::T, (0.0, Inf), "saturation loading")
    (K₀::T, (0.0, Inf), "affinity parameter")
    (E::T, (-Inf, 0.0), "energy parameter")
    (Bᵢᵩ::T, (-Inf, Inf), "adsorbate-adsorbent interaction coefficient")
end

function gibbs_excess_free_energy(model::ThermodynamicLangmuir, T::A, θ::AbstractVector{V}) where {V, A}
    _1 = one(eltype(T))
    Bᵢᵩ = model.Bᵢᵩ
    T⁻¹ = _1/T
    θᵢ, θᵩ = θ
    τᵢᵩ = Bᵢᵩ*T⁻¹
    Gᵢᵩ = exp(-0.3*τᵢᵩ)
    gᴱ_RT = θᵢ*θᵩ*τᵢᵩ*(Gᵢᵩ - _1)/(θᵢ*Gᵢᵩ + θᵩ)
    return gᴱ_RT
end

function thermodynamic_langmuir_activity_coefficient(model::ThermodynamicLangmuir, T, θᵢ)
    θᵩ = 1 - θᵢ
    Bᵢᵩ = model.Bᵢᵩ
    T⁻¹ = 1/T
    τᵢᵩ = Bᵢᵩ*T⁻¹
    Gᵢᵩ = exp(-0.3*τᵢᵩ)
    lnγᵢ = (Gᵢᵩ - 1)*τᵢᵩ*θᵩ*θᵩ/(θᵩ + θᵢ*Gᵢᵩ)^2
    lnγᵩ = Gᵢᵩ*(Gᵢᵩ - 1)*τᵢᵩ*θᵢ*θᵢ/(θᵩ + θᵢ*Gᵢᵩ)^2
    #lnγᵢ = θᵩ*θᵩ*τᵢᵩ*Gᵢᵩ*( lnγᵢ*(θᵢ + θᵩ*Gᵢᵩ)^-2 + (θᵩ + θᵢ*Gᵢᵩ)^-2)
    #lnγᵩ = θᵢ*θᵢ*τᵢᵩ*Gᵢᵩ*(Gᵢᵩ*(θᵩ + θᵢ*Gᵢᵩ)^-2 + (θᵢ + θᵩ*Gᵢᵩ)^-2)
    #exp((1 - Gᵢᵩ)*τᵢᵩ/(Gᵢᵩ)^2)
    return SVector((exp(lnγᵢ), exp(lnγᵩ)))
end

function activity_coefficient(model::ThermodynamicLangmuir, T::A, θ::AbstractVector{V}) where {V, A}
    fun = let model = model, T = T
    θ -> gibbs_excess_free_energy(model, T, θ)
    end
    #cfg = ForwardDiff.GradientConfig(fun, θ, autochunk(θ))
    cache = similar(θ)
    return exp.(ForwardDiff.gradient!(cache, fun, θ))
end

function loading_x0(model::ThermodynamicLangmuir, p, T)
    guess_model = from_vec(LangmuirS1, (model.M, model.K₀, model.E))
    return loading(guess_model, p, T)
end

function isotherm_res(model::ThermodynamicLangmuir, q, p, T::A) where {A}
    _1 = one(eltype(T))
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    θᵢ = q/M
    θᵩ = _1 - θᵢ
    γᵢ, γᵩ = thermodynamic_langmuir_activity_coefficient(model, T, θᵢ)
    return q - M * K *p / (γᵢ/γᵩ + K*p)
end

function henry_coefficient(model::ThermodynamicLangmuir, T)
    #=
    thermodynamic langmuir is a langmuir with Kx = K/dγ
    dγ = γᵢ/γᵩ(q)
    =#
    _0 = zero(eltype(T))
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    Bᵢᵩ = model.Bᵢᵩ
    T⁻¹ = 1/T
    τᵢᵩ = Bᵢᵩ*T⁻¹
    Gᵢᵩ = exp(-0.3*τᵢᵩ)
    dγ0 = exp((Gᵢᵩ - 1)*τᵢᵩ) #γᵢ/γᵩ when q = 0
    return M*K/dγ0

    #=
    q_0 = _0
    f = let model = model, T = T
        (∂q, ∂p) -> isotherm_res(model, ∂q, ∂p, T)
    end
    _f,_df = fgradf2(f, q_0, _0)
    ∂f0_∂q, ∂f0_∂P = _df
    return -∂f0_∂P/∂f0_∂q
    =#
end

function isosteric_heat(model::ThermodynamicLangmuir, p, T; Vᵃ = zero(eltype(p)), Vᵍ = Rgas(model)*T/p)
    q = loading(model, p, T)
    f = let model = model, q = q
        (∂p, ∂T) -> isotherm_res(model, q, ∂p, ∂T)
    end
    _f,_df = fgradf2(f, p, T)
    ∂n_∂p, ∂n_∂T = _df
    return T*(Vᵍ - Vᵃ)*∂n_∂T/∂n_∂p
end

function loading_impl(model::ThermodynamicLangmuir{L}, p, T) where L
    _1 = one(eltype(p))
    q0 = loading_x0(model, p, T)
    M = model.M
    K₀ = model.K₀
    E = model.E
    K = K₀*exp(-E/(Rgas(model)*T))
    f01(qq) = isotherm_res(model, qq, p, T)
    prob = Roots.ZeroProblem(to_halley(f01), q0)
    return Roots.solve(prob, Roots.Halley())
end

function loading(model::ThermodynamicLangmuir, p, T)
    return loading_impl(model, p, T)
end

requires_integration_sp_res(::ThermodynamicLangmuir) = Val{true}()

function pressure_x0(model::M, Π, T, ::typeof(sp_res)) where M <: ThermodynamicLangmuir
    γᵢ₀, γᵩ₀ = thermodynamic_langmuir_activity_coefficient(model, T, zero(Π))
    guess_model_0 = from_vec(LangmuirS1, SVector(promote(model.M, model.K₀*γᵩ₀/γᵢ₀, model.E)))
    p0 = pressure(guess_model_0, Π, T, sp_res)
    γᵢ∞, γᵩ∞ = thermodynamic_langmuir_activity_coefficient(model, T, oneunit(Π))
    guess_model_inf = from_vec(LangmuirS1, SVector(promote(model.M, model.K₀*γᵩ∞/γᵢ∞, model.E)))
    pinf = pressure(guess_model_inf, Π, T, sp_res)
    return 0.5*(p0 + pinf)
end

function saturated_loading(model::ThermodynamicLangmuir,T)
    return model.M*one(eltype(T))
end


function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: ThermodynamicLangmuir
    langmuir_model = x0_guess_fit(LangmuirS1,data)
    M, K₀, E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E    
    _0 = prevfloat(zero(M))
    return T(M, K₀, E, _0)
end

export ThermodynamicLangmuir, gibbs_excess_free_energy
