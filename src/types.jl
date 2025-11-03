"""
    FieldType

A type used to specify what type of physical field, such as traction or displacement.
"""
abstract type FieldType end

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

# Constructor without type parameters
function BoundaryData(field_type::F; 
        points = [zeros(3)], 
        interior_points = mean(points) |> x -> [SVector(x...)],
        fields = [p .* 0.0 for p in points]
    ) where F <: FieldType

    Dim = length(points[1])

    return BoundaryData{F,Dim}(interior_points, points, fields)
end