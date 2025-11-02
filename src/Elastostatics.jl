module Elastostatics

export DisplacementType,TractionType,BoundaryData, Elastic

using Accessors
using LinearAlgebra
using StaticArrays: SVector
using MultipleScattering
import Statistics: mean

# Write your package code here.

# for ploting recipes
using RecipesBase

include("types.jl")


end
