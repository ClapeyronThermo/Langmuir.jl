struct ExtendedLangmuir{T,𝕀} <: IsothermModel{T}
    isotherms::𝕀
end

struct IASTModels{T,𝕀} <: IsothermModel{T}
    isotherms::𝕀
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


function isosteric_heat(model::ExtendedLangmuir{_T, I}, pᵢ, T; Vᵃ = zeros(eltype(pᵢ)), Vᵍ = Rgas(model).*T./pᵢ) where {_T, I <: Tuple{Vararg{<:LangmuirS1{_T}}}}

    pᵢ_T = [pᵢ; T]
    f(pᵢ_T) = loading(model, first(pᵢ_T, length(pᵢ)), last(pᵢ_T))
    ȷ_I = ForwardDiff.jacobian(f, pᵢ_T)
    
    ∂pᵢ∂T = @views ȷ_I[:, 1:end - 1]\ȷ_I[:, end]
    return @. T*(Vᵍ - Vᵃ).*∂pᵢ∂T

end


function isosteric_heat(model::ExtendedLangmuir{_T, I}, p, T, yᵢ; Vᵃ = zeros(eltype(p.*yᵢ)), Vᵍ = Rgas(model).*T./(p.*yᵢ)) where {_T, I <: Tuple{Vararg{<:LangmuirS1{_T}}}}
    _y = yᵢ/sum(yᵢ) 
    pᵢ = p*_y
    isosteric_heat(model, pᵢ, T, Vᵃ = Vᵃ, Vᵍ = Vᵍ)
end




export ExtendedLangmuir, IASTModels