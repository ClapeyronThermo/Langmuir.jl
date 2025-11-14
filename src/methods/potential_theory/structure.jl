struct ExternalField1DCartesian{I <: Int, F <: AbstractPotential} <: PTAStructure1D
    potential::F
    ngrids::Tuple{I}
end

function upper_bound(potential::DRA)
    z0 = potential.z0
    return z0 
end

function upper_bound(potential::MultiComponentDRA)
    potentials = potential.potentials
    z0_values = map(p -> upper_bound(p), potentials)
    return maximum(z0_values)
end

function lower_bound(potential::DRA; cutoff = 1e-8)
    return potential.z0*cutoff
end

function lower_bound(potential::MultiComponentDRA; cutoff = 1e-8)
    potentials = potential.potentials
    z0_values = map(p -> lower_bound(p; cutoff = cutoff), potentials)
    return minimum(z0_values)
end

function bounds(potential::DRA)
    return [lower_bound(potential), upper_bound(potential)]
end

function ExternalField1DCartesian(potential::F, ngrids = (101, )) where {F <: AbstractPotential}
    return ExternalField1DCartesian(potential, ngrids)
end

function create_grid(structure::S) where {S <: PTAStructure1D}
    potential = structure.potential
    return LinRange(lower_bound(potential),
                    upper_bound(potential),
                    structure.ngrids[1])
end

export ExternalField1DCartesian