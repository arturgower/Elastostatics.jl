module MethodOfFundamentalSolutions

export PointCloud, outward_normals, points_in_shape
export FieldResults, field
export DisplacementType, TractionType, BoundaryData, Elastic, Acoustic

using Accessors
using LinearAlgebra
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