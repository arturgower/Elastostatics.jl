module MethodOfFundamentalSolutions

using Accessors
using LinearAlgebra
using BlockArrays: mortar
using StaticArrays: SVector
using MultipleScattering

import MultipleScattering: PhysicalMedium, ScalarMedium, spatial_dimension, field_dimension, Acoustic, Shape, Box, bounding_box, points_in_shape, cartesian_to_radial_coordinates, radial_to_cartesian_transform, cartesian_to_radial_transform, field
export cartesian_to_radial_coordinates, radial_to_cartesian_transform, cartesian_to_radial_transform

import Statistics: mean

# for ploting recipes
using RecipesBase


export BoundaryData, outward_normals, points_in_shape
include("boundarydata.jl")

export FieldResult, greens, field, source_positions
export source_system, FundamentalSolution
include("solve.jl")

export DisplacementType, TractionType, Elastostatic, Acoustic
include("physics/elastic.jl")
include("../plot/plot.jl")

end