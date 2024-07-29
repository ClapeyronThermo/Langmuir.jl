struct HenryIsotherm{T} <: IsothermModel{T}
    Kh::T
end

function sp_res(model::HenryIsotherm,p)
    Kh = model.Kh
    return Kh*p
end

function henry_coefficient(model::HenryIsotherm)
    return model.Kh
end

sp_res_pressure_impl(model::HenryIsotherm,Π) = Π/model.Kh

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: HenryIsotherm
    # use first two data points to get the slope
    l,p = data.l,data.p
    return HenryIsotherm(p\l)
end

export HenryIsotherm
