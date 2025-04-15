abstract type MultiComponentIsothermModel{T} end

struct ExtendedLangmuir{T,𝕀} <: MultiComponentIsothermModel{T}
    isotherms::𝕀
end

struct IASTModels{T,𝕀} <: MultiComponentIsothermModel{T}
    isotherms::𝕀
end

struct ThermodynamicLangmuirModels{T,𝕀} <: MultiComponentIsothermModel{T}
    isotherms::𝕀
end

struct aNRTLModel{𝕀, T <: Array}
    isotherms::𝕀
    Βᵢⱼ::T
end

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

function _tlangmuirmodels(::Type{T}, isotherms::𝕀) where {T,𝕀}
    return ThermodynamicLangmuirModels{T,𝕀}(isotherms)
end

_tlangmuirmodels(isotherms::I) where I = _tlangmuirmodels(eltype(first(isotherms)), isotherms)

function ThermodynamicLangmuirModels(m_first::I, m_rest::Vararg{I}) where {I <: ThermodynamicLangmuir}
    return _tlangmuirmodels((m_first, m_rest...))
end

function aNRTLModel(I::N) where N <: ThermodynamicLangmuirModels
    isotherms = I.isotherms
    Bᵢᵩ = getfield.(isotherms, :Bᵢᵩ) |> collect
    Bᵢⱼ = Bᵢᵩ .- Bᵢᵩ' # Estimated interaction parameters from pure isotherm 
    return aNRTLModel(isotherms, Bᵢⱼ)
end

function aNRTLModel(I::N) where N <: MultiComponentIsothermModel
    isotherms = I.isotherms
    NIsotherms = length(isotherms)
    Bᵢᵩ = zeros(NIsotherms)
    Bᵢⱼ = Bᵢᵩ .- Bᵢᵩ' # Estimated interaction parameters from pure isotherm 
    return aNRTLModel(I, Bᵢⱼ)
end

function loading(model::ExtendedLangmuir{_T, I}, p, T, yᵢ) where {_T, I <: Tuple{Vararg{<:LangmuirS1{_T}}}}
    _y = yᵢ/sum(yᵢ) 
    pᵢ = p*_y
    return loading(model, pᵢ, T)
end

function loading(model::ExtendedLangmuir{_T, I}, pᵢ, T) where {_T, I <: Tuple{Vararg{<:LangmuirS1{_T}}}}

    _1_∑kP = one(eltype(model))
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


function loading(model::ExtendedLangmuir{_T, I}, p, T, y) where {_T, I <: Tuple{Vararg{<:MultiSite{_T}}}}
        
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


function isosteric_heat(model::ExtendedLangmuir{_T, I}, pᵢ, T) where {_T, I <: Tuple{Vararg{<:LangmuirS1{_T}}}}

#=     pᵢ_T = [pᵢ; T]
    f(pᵢ_T) = loading(model, first(pᵢ_T, length(pᵢ)), last(pᵢ_T))
    ȷ_I = ForwardDiff.jacobian(f, pᵢ_T)
    
    ∂pᵢ∂T = @views ȷ_I[:, 1:end - 1]\ȷ_I[:, end]
    return @. T*(Vᵍ - Vᵃ).*∂pᵢ∂T #No excess heat due to mixing. =#

    return map(p_model -> isosteric_heat(last(p_model), first(p_model), T), zip(pᵢ, model.isotherms))

end

function isosteric_heat(model::ExtendedLangmuir{_T, I}, pᵢ, T) where {_T, I <: Tuple{Vararg{<: MultiSite{_T} }}}
    
        return map(p_model -> isosteric_heat(last(p_model), first(p_model), T), zip(pᵢ, model.isotherms))
    
end


function isosteric_heat(model::ExtendedLangmuir{_T, I}, p, T, yᵢ) where {_T, I <: Tuple{Vararg{ <: Union{MultiSite{_T}, LangmuirS1{_T}}}}}
    _y = yᵢ/sum(yᵢ) 
    pᵢ = p*_y
    isosteric_heat(model, pᵢ, T)
end


function gibbs_excess_free_energy(model::aNRTLModel, T, x::Vector)
    
    θ = x./sum(x)
    Bᵢⱼ = model.Βᵢⱼ
    T⁻¹ = one(T)./T

    τᵢⱼ = Bᵢⱼ*T⁻¹
    Gᵢⱼ = exp.(-0.3.*τᵢⱼ)

    ∑ⱼθⱼτᵢⱼ = τᵢⱼ*θ
    ∑ₖθₖGₖᵢ = Gᵢⱼ'*θ

    gᴱ_RT = sum((∑ⱼθⱼτᵢⱼ./∑ₖθₖGₖᵢ).*θ)

    return gᴱ_RT
    
end

function activity_coefficient(model::aNRTLModel, T, x::Vector)
    
    fun(x) = gibbs_excess_free_energy(model, T, x)

    return exp.(gradient(fun, x))
end

export ExtendedLangmuir, IASTModels, ThermodynamicLangmuirModels, aNRTLModel