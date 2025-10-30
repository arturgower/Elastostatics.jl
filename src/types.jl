## field types
"""
    FieldType

A type used to specify what type of physical field, such as traction or displacement.
"""
abstract type FieldType end

struct DisplacementType <: FieldType end
struct TractionType <: FieldType end

struct BoundaryData{F <: FieldType,Dim}
    points::Vector{SVector{Dim,Float64}}
    fields::Vector{SVector{Dim,Float64}}
end

import ElasticWaves: Elastic