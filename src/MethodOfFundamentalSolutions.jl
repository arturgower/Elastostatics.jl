module MethodOfFundamentalSolutions

export BoundaryData, outward_normals, points_in_shape
export FieldResults, greens, field, source_positions
export DisplacementType, TractionType, Elastostatic, Acoustic

using Accessors
using LinearAlgebra
using BlockArrays: mortar
using StaticArrays: SVector
using MultipleScattering

import MultipleScattering: PhysicalMedium, ScalarMedium, Acoustic, Shape, Box, bounding_box, points_in_shape, field

import Statistics: mean

# for ploting recipes
using RecipesBase

include("boundarydata.jl")
include("solve.jl")

include("physics/elastic.jl")
include("../plot/plot.jl")

end