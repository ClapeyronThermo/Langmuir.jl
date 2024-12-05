module Langmuir

using NLSolvers
using LinearAlgebra
using ForwardDiff
using LogExpFunctions
using StaticArrays
using Integrals
using Tables, TableOperations
using Roots: Roots, solve
import PolyLog
import CommonSolve
using BlackBoxOptim

abstract type IsothermModel{T} end

include("utils.jl")
include("isotherm_data.jl")
include("models/models.jl")
include("methods/methods.jl")

end
