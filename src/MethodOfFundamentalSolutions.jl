module MethodOfFundamentalSolutions

using Accessors
using LinearAlgebra
using BlockArrays: mortar
using StaticArrays: SVector
using MultipleScattering

import MultipleScattering: PhysicalMedium, ScalarMedium, spatial_dimension, field_dimension, Acoustic, Shape, Box, bounding_box, points_in_shape, cartesian_to_radial_coordinates, field

import Statistics: mean

# for ploting recipes
using RecipesBase


export BoundaryData, outward_normals, points_in_shape
include("boundarydata.jl")

export FieldResults, greens, field, source_positions
export FundamentalSolution
include("solve.jl")

export DisplacementType, TractionType, Elastostatic, Acoustic
include("physics/elastic.jl")
include("../plot/plot.jl")

end