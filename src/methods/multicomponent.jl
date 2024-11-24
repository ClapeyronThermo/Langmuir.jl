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


function loading(model::ExtendedLangmuir{_T, I}, p, T, y) where {_T, I <: Tuple{Vararg{<:LangmuirS1{_T}}}}

    _y = y/sum(y) 
    p_i = p*_y
    one_plus_Kp = one(eltype(model))
    loadings = similar(y, eltype(model))
    models = model.isotherms

    for i in eachindex(y)
            loadings[i] = begin
                M, K₀, E = models[i].M, models[i].K₀, models[i].E
                MKpy = M*K₀*exp(-E/(Rgas(model)*T))*p_i[i]
                one_plus_Kp += MKpy/M
                MKpy
            end
    end

    return loadings./one_plus_Kp
end

function loading(model::IASTModels, p, T, y; method = IASTNestedLoop(), gas_model = nothing, x0 = nothing, maxiters = 100, reltol = 1e-12, abstol = 1e-10)
     
    (n_total, x, is_success) = iast(model.isotherms, p, T, y, method, gas_model, x0 = x0, maxiters = maxiters, reltol = reltol, abstol = abstol)

    if is_success == :success
        loadings = x.*n_total
    else
        # Convergence failed
        error("Convergence failed - current number of iterations is $maxiters, consider increasing to meet tolerances.")
    end

    return loadings
end

export ExtendedLangmuir, IASTModels