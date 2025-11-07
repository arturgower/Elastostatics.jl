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
export particular_solution
include("boundarydata.jl")

export FieldResult, FundamentalSolution
export greens, field, source_positions
export solve, TikhonovSolver, system_matrix
include("solve.jl")

export DisplacementType, TractionType, Elastostatic, Acoustic
export ParticularGravity
include("physics/elastic.jl")
include("../plot/plot.jl")

end