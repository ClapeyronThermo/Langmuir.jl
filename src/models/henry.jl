struct Henry{T} <: IsothermModel{T}
    Kh::T
end

function sp_res(model::Henry,p)
    Kh = model.Kh
    return Kh*p
end

function henry_coefficient(model::Henry)
    return model.Kh
end

pressure_impl(model::Henry,Π,T,::typeof(sp_res),approx) = Π/model.Kh

function x0_guess_fit(::Type{T},data::AdsIsoTData) where T <: Henry
    # use first two data points to get the slope
    l,p = data.l,data.p
    return Henry(p\l)
end

export Henry
