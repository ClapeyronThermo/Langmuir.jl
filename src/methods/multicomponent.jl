abstract type MultiComponentIsothermModel{T} end

struct ExtendedLangmuir{T,𝕀} <: MultiComponentIsothermModel{T}
    isotherms::𝕀
end

struct IASTModels{T,𝕀} <: MultiComponentIsothermModel{T}
    isotherms::𝕀
end

struct aNRTLModel{T, 𝕀, B <: AbstractMatrix{T}} <: MultiComponentIsothermModel{T}
    isotherms::𝕀
    Βᵢⱼ::B
end

const LangmuirS1Tuple{T,N} = NTuple{N,LangmuirS1{T}} where {T,N}
const MultiSiteTuple{T,N,I} = NTuple{N,MultiSite{T,I}} where {T,N,I}

function _extendedlangmuir(::Type{T}, isotherms::𝕀) where {T,𝕀}
    return ExtendedLangmuir{T,𝕀}(isotherms)
end

_extendedlangmuir(isotherms::I) where I = _extendedlangmuir(eltype(first(isotherms)), isotherms)

function ExtendedLangmuir(m_first::I, m_rest::Vararg{I}) where {I <: Union{LangmuirS1, MultiSite}}
    return _extendedlangmuir((m_first, m_rest...))
end


function _iastmodels(::Type{T}, isotherms::𝕀) where {T,𝕀}
    return IASTModels{T,𝕀}(isotherms)
end

_iastmodels(isotherms::I) where I = _iastmodels(eltype(first(isotherms)), isotherms)

function IASTModels(m_first::I, m_rest::Vararg{I}) where {I <: IsothermModel}
    return _iastmodels((m_first, m_rest...))
end

function aNRTLModel(models::Tuple{I, Vararg{I}}) where I <: ThermodynamicLangmuir
    Bᵢᵩ = getfield.(models, :Bᵢᵩ) |> collect
    Bᵢⱼ = Bᵢᵩ .- Bᵢᵩ' # Estimated interaction parameters from pure isotherm
    T = eltype(Bᵢⱼ)
    return aNRTLModel{T,typeof(models),typeof(Bᵢⱼ)}(models, Bᵢⱼ)
end

function aNRTLModel(models::Tuple{I, Vararg{I}}) where I <: IsothermModel
    isotherms = collect(models)
    NIsotherms = length(isotherms)
    T = Base.promote_eltype(models...)
    Bᵢᵩ = zeros(T,(NIsotherms,NIsotherms)) #Ideal model
    return aNRTLModel{T,typeof(models),typeof(Bᵢⱼ)}(models, Bᵢⱼ)
end

function loading(model::ExtendedLangmuir{_T, I}, p, T, yᵢ) where {_T, I <: LangmuirS1Tuple{_T}}
    _y = yᵢ/sum(yᵢ)
    pᵢ = p*_y
    return loading(model, pᵢ, T)
end

function loading(model::ExtendedLangmuir{_T, I}, pᵢ, T) where {_T, I <: LangmuirS1Tuple{_T}}

    _1_∑kP = one(eltype(T))
    loadings = similar(pᵢ)
    models = model.isotherms

    for i ∈ eachindex(pᵢ)
            loadings[i] = begin
                M, K₀, E = models[i].M, models[i].K₀, models[i].E
                MKpy = M*K₀*exp(-E/(Rgas(model)*T))*pᵢ[i]
                _1_∑kP += MKpy/M
                MKpy
            end
    end

    return loadings./_1_∑kP
end

function loading(model::ExtendedLangmuir{_T, I}, p, T, y) where {_T, I <: MultiSiteTuple{_T}}

        unpack_multisite = getfield.(model.isotherms, :isotherms)

        Es = collect.(map(x -> getfield.(x, :E), unpack_multisite) |> collect)

        N_sites =  length(Es[1])

        sites_to_match = sortperm.(Es)

        models = [ExtendedLangmuir(
            getindex.(collect(unpack_multisite), getindex.(sites_to_match, i))...
            )
             for i in 1:N_sites]

        loadings = mapreduce(m -> loading(m, p, T, y), +, models)

        return loadings
end


function loading(model::IASTModels, p, T, y; method = IASTNestedLoop(), gas_model = nothing, x0 = nothing, maxiters = 100, reltol = 1e-12, abstol = 1e-10)

    (n_total, x, is_success) = iast(model.isotherms, p, T, y, method, gas_model, x0 = x0, maxiters = maxiters, reltol = reltol, abstol = abstol)

    loadings = similar(y)

    if is_success == :success
        loadings .= x.*n_total
    else
        # Convergence failed
        error("Convergence failed - current number of iterations is $maxiters, consider increasing to meet tolerances.")
    end

    return loadings
end


function isosteric_heat(model::IASTModels, p, T, y; Vg = Rgas(model).*T./(p.*y), Vₐ = zeros(length(y)), method = IASTNestedLoop(), gas_model = nothing, x0 = nothing, maxiters = 100, reltol = 1e-12, abstol = 1e-10)

    (n_total, x, is_success) = iast(model.isotherms, p, T, y, method, gas_model, x0 = x0, maxiters = maxiters, reltol = reltol, abstol = abstol)

    if is_success == :success

        Pᵢ₀_π = y.*p./x
    else
        # Convergence failed
        error("Convergence failed - current number of iterations is $maxiters, consider increasing to meet tolerances.")
    end

    ∂ψ_∂p = map(m_Pᵢ₀_π -> -1.0 * loading(first(m_Pᵢ₀_π), last(m_Pᵢ₀_π), T)/last(m_Pᵢ₀_π), zip(model.isotherms, Pᵢ₀_π))

    ∂ψ_∂T = map(m_Pᵢ₀_π -> -1.0 * ForwardDiff.derivative(T -> sp_res(first(m_Pᵢ₀_π), last(m_Pᵢ₀_π), T),  T), zip(model.isotherms, Pᵢ₀_π))

    return T.*(Vg .- Vₐ).*(∂ψ_∂T./∂ψ_∂p).*x
end


function isosteric_heat(model::ExtendedLangmuir{_T, I}, pᵢ, T) where {_T, I <: LangmuirS1Tuple{_T}}

#=     pᵢ_T = [pᵢ; T]
    f(pᵢ_T) = loading(model, first(pᵢ_T, length(pᵢ)), last(pᵢ_T))
    ȷ_I = ForwardDiff.jacobian(f, pᵢ_T)

    ∂pᵢ∂T = @views ȷ_I[:, 1:end - 1]\ȷ_I[:, end]
    return @. T*(Vᵍ - Vᵃ).*∂pᵢ∂T #No excess heat due to mixing. =#

    return map(p_model -> isosteric_heat(last(p_model), first(p_model), T), zip(pᵢ, model.isotherms))

end

function isosteric_heat(model::ExtendedLangmuir{_T, I}, pᵢ, T) where {_T, I <: MultiSiteTuple{_T}}
    return map(p_model -> isosteric_heat(last(p_model), first(p_model), T), zip(pᵢ, model.isotherms))
end


function isosteric_heat(model::ExtendedLangmuir{_T, I}, p, T, yᵢ) where {_T, I <: Union{LangmuirS1Tuple{_T},MultiSiteTuple{_T}}}
    pᵢ = p*_y
    isosteric_heat(model, pᵢ, T)
end

function gibbs_excess_free_energy(model::aNRTLModel, T, x::AbstractVector)
    Bᵢⱼ = model.Βᵢⱼ
    i1,i2 = size(Bᵢⱼ)
    _0 = zero(Base.promote_eltype(Bᵢⱼ,T,x))
    ∑x⁻¹ = 1/sum(x)
    T⁻¹ = one(T)/T
    gᴱ_RT = _0
    for i ∈ 1:i1
        ∑τⱼᵢGⱼᵢxⱼ = _0
        ∑Gⱼᵢxⱼ = _0
        xᵢ = x[i]*∑x⁻¹
        for j ∈ 1:i2
            xⱼ = x[j]*∑x⁻¹
            τⱼᵢ = Bᵢⱼ[j,i]*T⁻¹
            Gⱼᵢ = exp(-0.3*τⱼᵢ)
            Gⱼᵢxⱼ = xⱼ*Gⱼᵢ
            ∑Gⱼᵢxⱼ += Gⱼᵢxⱼ
            ∑τⱼᵢGⱼᵢxⱼ += Gⱼᵢxⱼ*τⱼᵢ
        end
        gᴱ_RT += xᵢ*∑τⱼᵢGⱼᵢxⱼ/∑Gⱼᵢxⱼ
    end
    return gᴱ_RT
#=
    θ = x./sum(x)
    Bᵢⱼ = model.Βᵢⱼ
    T⁻¹ = one(T)./T

    τᵢⱼ = Bᵢⱼ*T⁻¹
    Gᵢⱼ = exp.(-0.3.*τᵢⱼ)

    ∑ⱼθⱼτᵢⱼ = τᵢⱼ*θ
    ∑ₖθₖGₖᵢ = Gᵢⱼ'*θ

    gᴱ_RT = sum((∑ⱼθⱼτᵢⱼ./∑ₖθₖGₖᵢ).*θ)

    return gᴱ_RT
    =#
end

function activity_coefficient(model::MultiComponentIsothermModel, T, x::AbstractVector)
    
    fun = let model = model, T = T
        x -> gibbs_excess_free_energy(model, T, x)
    end
    cache = similar(x)
    return exp.(ForwardDiff.gradient!(cache, fun, x))
end


export ExtendedLangmuir, IASTModels, ThermodynamicLangmuirModels, aNRTLModel