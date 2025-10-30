"""
    FieldType

A type used to specify what type of physical field, such as traction or displacement.
"""
abstract type FieldType end

struct DisplacementType <: FieldType end
struct TractionType <: FieldType end

"""
BoundaryData{F,Dim}

Type parameters
- F <: FieldType: the physical field witehr displacement or traction.
- Dim: Integer-sized compile-time dimension of the spatial dimension

Fields
- `interior_points::Vector{SVector{Dim,Float64}}`
    A vector of spatial points in the interior of the domain used to determine what is inside the domain. See function.... 
- `points::Vector{SVector{Dim,Float64}}`
    A vector of spatial points on the boundary
- `fields::Vector{SVector{Dim,Float64}}`
    `fields[i]` is the value of the physical field at `points[i]`

Notes
- It is expected that `length(points) == length(fields)` and that entries are aligned by index.
"""
struct BoundaryData{F <: FieldType,Dim}
    interior_points::Vector{SVector{Dim,Float64}}
    points::Vector{SVector{Dim,Float64}}
    fields::Vector{SVector{Dim,Float64}}
end

import ElasticWaves: Elastic