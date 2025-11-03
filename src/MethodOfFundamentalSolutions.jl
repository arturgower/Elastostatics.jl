module MethodOfFundamentalSolutions

export DisplacementType,TractionType,BoundaryData, Elastic, Acoustic

using Accessors
using LinearAlgebra
using StaticArrays: SVector
using MultipleScattering

import MultipleScattering: PhysicalMedium, ScalarMedium, Acoustic, Shape, Box, bounding_box
import Statistics: mean

# for ploting recipes
using RecipesBase

include("types.jl")
include("pointcloud.jl")
include("physics/elastic.jl")


end
