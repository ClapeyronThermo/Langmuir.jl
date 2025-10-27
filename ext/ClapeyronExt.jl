module ClapeyronExt

using Langmuir, Clapeyron

function Langmuir.chemical_potential(model::M, p, T, z, phase = :stable) where {M <: Clapeyron.EoSModel}
    return Clapeyron.chemical_potential(model, p, T, z, phase = phase)
end

function Langmuir.VT_chemical_potential(model::M, v, T, z::Z) where {M <: Clapeyron.EoSModel, Z <: Number}
    if length(z) == 1
        z = [z]
    end
    return first(Clapeyron.VT_chemical_potential(model, v, T, z))
end

Langmuir.Rg(eos::M) where M <: Clapeyron.EoSModel = 8.31446261815324 

function Langmuir.VT_chemical_potential(model::M, v, T, z::Z) where {M <: Clapeyron.EoSModel, Z <: AbstractVector}
    return Clapeyron.VT_chemical_potential(model, v, T, z)
end

function Langmuir.idealmodel(model::M) where {M <: Clapeyron.EoSModel}
    return Clapeyron.idealmodel(model)
end

function Langmuir.molar_density(model::M, p, T, z) where {M <: Clapeyron.IdealModel}
    return Clapeyron.molar_density(model, p, T, z)
end

function Langmuir.volume(model::M, P, T, z, phase = :stable) where {M <: Clapeyron.EoSModel}
    return Clapeyron.volume(model, P, T, z, phase = phase)
end

function Langmuir.VT_isstable(model::M, v, T, z::Z) where {M <: Clapeyron.EoSModel, Z <: Number}
    z = [z]
    return Clapeyron.VT_isstable(model, v, T, z)
end

function Langmuir.VT_isstable(model::M, v, T, z::Z) where {M <: Clapeyron.EoSModel, Z <: AbstractVector}
    return Clapeyron.VT_isstable(model, v, T, z)
end

end #module