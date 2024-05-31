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