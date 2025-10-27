abstract type AbstractPTAStructure end
abstract type PTAStructure1D <: AbstractPTAStructure end
abstract type AbstractPotential end
abstract type SingleComponentPotential <: AbstractPotential end
abstract type MultiComponentPotential <: AbstractPotential end

struct PTASystem{M, S <: AbstractPTAStructure} 
    model::M
    structure::S
    function PTASystem(model::M, structure::S) where {M, S <: AbstractPTAStructure}
        length(model) == length(structure.potential) ||
            throw(ArgumentError("model and structure.potential must have the same length"))
    new{typeof(model),S}(model, structure)
    end
end

function PTASystem(model, potential::P) where {P <: AbstractPotential}
    structure = ExternalField1DCartesian(potential)
    return PTASystem(model, structure)
end

function Rg(model::T) where T 
    return 8.31446261815324  # Universal gas constant in J/(mol·K)
end


export PTASystem