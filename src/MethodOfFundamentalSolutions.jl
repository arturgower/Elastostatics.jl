module MethodOfFundamentalSolutions

export PointCloud, outward_normals, outward_normals_and_neighbour_distance, points_in_shape
export FieldResults, greens, field, source_positions
export DisplacementType, TractionType, BoundaryData, Elastic, Acoustic

using Accessors
using LinearAlgebra
using BlockArrays: mortar
using StaticArrays: SVector
using MultipleScattering

import MultipleScattering: PhysicalMedium, ScalarMedium, Acoustic, Shape, Box, bounding_box, points_in_shape, field

import Statistics: mean

# for ploting recipes
using RecipesBase

include("types.jl")
include("pointcloud.jl")
include("physics/elastic.jl")
include("../plot/plot.jl")

end