module AdsorbedSolutionTheory

using NLSolvers
using LinearAlgebra
using ForwardDiff
using LogExpFunctions
using StaticArrays
using Tables
using Roots

include("utils.jl")
include("isotherm_data.jl")
include("methods/methods.jl")
include("models/models.jl")
# Write your package code here.

end
