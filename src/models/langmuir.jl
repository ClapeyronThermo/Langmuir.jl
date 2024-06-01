#general multisite langmuir.
struct Langmuir{T} <: IsothermModel{T}
    M::T
    K::T
end

function sp_res(model::Langmuir,p)
    M = model.M
    K = model.K
    return M*log1p(K*p)
end

#optimizations for Langmuir, not necessary, but improve performance
henry_coefficient(model::Langmuir) = model.M*model.K
saturated_loading(model::Langmuir) = model.M
sp_res_pressure_impl(model::Langmuir,Π) = expm1(Π/model.M)/model.K

struct DualSiteLangmuir{T} <: IsothermModel{T}
    M1::T
    K1::T
    M2::T
    K2::T
end

function sp_res(model::DualSiteLangmuir,p)
    M1,K1,M2,K2 = model.M1,model.K1,model.M2,model.K2
    return M1*log1p(K1*p) + M2*log1p(K2*p)
end

#optimizations for DualSiteLangmuir, not necessary, but improve performance
henry_coefficient(model::DualSiteLangmuir) = model.M1*model.K1 + model.M2*model.K2
saturated_loading(model::DualSiteLangmuir) = model.M1 + model.M2

export Langmuir, DualSiteLangmuir