#useful for defining carrier gas

struct ZeroIsotherm{T} <: IsothermModel{T} end


"""
    `ZeroIsotherm()`

    ZeroIsotherm <: IsothermModel

ZeroIsotherm represents an isotherm that has no adsorption. it can be used to define carrier gases.

## Description

A zero isotherm is defined by:

n == 0

where:
- n is the loading of the adsorbate on the adsorbent
"""
ZeroIsotherm() = ZeroIsotherm{Float64}()
ZeroIsotherm(T) = ZeroIsotherm{T}()

sp_res(model::ZeroIsotherm, p, T) = zero(Base.promote_eltype(model,p,T))
loading(model::ZeroIsotherm, p, T) = zero(Base.promote_eltype(model,p,T))
henry_coefficient(model::ZeroIsotherm, T) = zero(Base.promote_eltype(model,T))
Base.iszero(::ZeroIsotherm) = true
