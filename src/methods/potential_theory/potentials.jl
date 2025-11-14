struct DRA{T} <: SingleComponentPotential
    ε0::T
    z0::T
    Β::T
end

struct MultiComponentDRA{T, I} <: MultiComponentPotential
    potentials::I
end

_multicomponentdra(::Type{T}, potentials::I) where {T, I} = MultiComponentDRA{T, I}(potentials)

_multicomponentdra(potentials::I) where I = _multicomponentdra(eltype(first(potentials)), potentials)

function MultiComponentDRA(m_first::I, m_rest::Vararg{I}) where I <: DRA
    return _multicomponentdra((m_first, m_rest...))
end

Base.length(model::DRA) = 1

Base.length(model::MultiComponentDRA) = length(model.potentials)

function potential(model::DRA, z::T) where T <: Real
    _1 = one(eltype(model))
    ε0 = model.ε0
    z0 = model.z0
    Β = model.Β
    _1_Β = _1 / Β
    _0 = zero(eltype(model))
    return ifelse(z < z0, ε0*abs(log(z0/z))^(_1_Β), _0)
end

function potential(model::MultiComponentDRA, z::T) where T <: Real
    
    f = let z = z
        p -> potential(p, z)
    end
    
    return [map(f, model.potentials)...]
end

export DRA, MultiComponentDRA, potential 

Base.eltype(potential::DRA{T}) where T = T
Base.eltype(potential::MultiComponentDRA{T, I}) where {T, I} = T