module AdsorbedSolutionTheory

using NLSolvers
using LinearAlgebra
using ForwardDiff
using LogExpFunctions
using StaticArrays
using Integrals
using Tables
using Roots: Roots, solve
import PolyLog
import CommonSolve

abstract type IsothermModel{T} end

include("utils.jl")
include("isotherm_data.jl")
include("methods/methods.jl")
include("models/models.jl")

end
