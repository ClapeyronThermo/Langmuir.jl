"""
    Quadratic(K₀a, K₀b, M, Ea, Eb)

    Quadratic <: IsothermModel

## Inputs

- `K₀a::T`: Affinity parameter A, `[1/Pa]`
- `K₀b::T`: Affinity parameter B, `[1/Pa^2]`
- `M::T`: Saturation loading, `[mol/kg]`
- `Ea::T`: Adsorption energy A, `[J/mol]`
- `Eb::T`: Adsorption energy B, `[J/mol]`

## Description

The `Quadratic` isotherm model is given by:

n = M × (Ka + 2Kb × p) × p / (1 + p × (Ka + Kb × p))

The model assumes that the affinity parameters `Ka` and `Kb` are temperature-dependent and follow the relation:

Ka = K₀a * exp(-Ea / (RT))

Kb = K₀b * exp(-Eb / (RT))

Where:
- `Ka` and `Kb` are the affinity parameters at temperature `T`,
- `R` is the gas constant,
- `T` is the absolute temperature.

"""
@with_metadata struct Quadratic{T} <: IsothermModel{T}
    (K₀a::T, (0.0, Inf), "Affinity parameter A")
    (K₀b::T, (0.0, Inf), "Affinity parameter B")
    (M::T, (0.0, Inf), "Saturation loading")
    (Ea::T, (-Inf, 0.0), "Energy parameter A")
    (Eb::T, (-Inf, 0.0), "Energy parameter B")
end

function sp_res(model::Quadratic, p, T)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    return M*log1p(p*(Ka + Kb*p))
end

function loading(model::Quadratic, p, T)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    _1 = one(eltype(p))
    return M*(Ka + 2.0*Kb*p)*p/(_1 + p*(Ka + Kb*p))
end

function pressure_impl(model::Quadratic, Π, T, ::typeof(sp_res), approx)
    K₀a, K₀b, M, Ea, Eb = model.K₀a, model.K₀b, model.M, model.Ea, model.Eb
    Ka = K₀a*exp(-Ea/(Rgas(model)*T))
    Kb = K₀b*exp(-Eb/(Rgas(model)*T))
    Kab = Ka/Kb
    return -0.5*Kab + sqrt(0.25*Kab*Kab + expm1(Π/M)/Kb)
end

function x0_guess_fit(::Type{T}, data::AdsIsoTData) where T <: Quadratic
    #langmuir_model = x0_guess_fit(Langmuir,data)
    #M, K₀, E = langmuir_model.M, langmuir_model.K₀, langmuir_model.E

    #l*(1 + p*(Ka + Kb*p)) = M*(Ka*p + 2*Kb*p*p)
    #p*Ka*l + Kb*p*p*l - M*Ka*p - 2*M*Kb*p*p = -l
    #-p*Ka*l - Kb*p*p*l + M*Ka*p + 2*M*Kb*p*p = l
    #-Ka*(p*l) -Kb*(p*p*l) + M*Ka*(p) + 2*M*Kb*(p*p) = l

    # Split data by temperature
    Ts, l_p = split_data_by_temperature(data)

    # Initialize vectors for Ka, Kb, and M values
    Kas = Vector{eltype(Ts)}(undef, length(l_p))
    Kbs = Vector{eltype(Ts)}(undef, length(l_p))
    Ms = Vector{eltype(Ts)}(undef, length(l_p))
   
    # Perform fitting for each (l, p) tuple
    for i in 1:length(l_p)
        l_i, p_i = l_p[i]
        Kaneg, Kbneg, MKa, MKb2 = hcat(p_i .* l_i, p_i .* p_i .* l_i, p_i, p_i .* p_i) \ l_i

        Ka = abs(Kaneg)
        Kb = abs(Kbneg)
        Ma = MKa / Ka
        Mb = abs(0.5 * MKb2 / Kb)
        Ms[i] = 0.5 * (Ma + Mb)
        Kas[i] = Ka
        Kbs[i] = Kb
    end

    M = sum(Ms)/length(Ms) #Mean of all values

    _1 = one(eltype(Ts))
    _1s = ones(eltype(Ts), length(Ts))

    if length(l_p) > 1
        logKa, Ea = hcat(_1s, _1s./ (Rgas(T).*Ts)) \ log.(Kas)
        Ka = exp(logKa)
        logKb, Eb = hcat(_1s, _1s./(Rgas(T).*Ts)) \ log.(Kbs)
        Kb = exp(logKb)
    else
        Ka = first(Kas)
        Kb = first(Kbs)
        Ea = _1
        Eb = _1
    end

    Quadratic(Ka, Kb, M, -Ea, -Eb)
end

export Quadratic
