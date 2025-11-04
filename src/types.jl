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

"""
    FieldResults{T<:Real,Dim,FieldDim}

Struct to hold results of simulations evaluated over different spatial positions.

# Type Parameters
- `T`: The numeric type for spatial coordinates (e.g., Float64)
- `Dim`: Spatial dimension of the problem (typically 2 or 3)
- `FieldDim`: Dimension of the field values (e.g., 1 for scalar fields, Dim for vector fields)

# Fields
- `x::Vector{SVector{Dim,T}}`: Vector of spatial positions
- `field::Vector{SVector{FieldDim,Complex{T}}}`: Vector of field values at corresponding positions

# Examples
```julia
# Create field results for a 2D problem with complex scalar field
positions = [SVector(0.0, 0.0), SVector(1.0, 0.0)]
fields = [SVector(1.0 + 0.0im), SVector(0.0 + 1.0im)]
results = FieldResults(positions, fields)
```
"""
struct FieldResults{T<:Real,Dim,FieldDim}
    "Spatial positions where field is evaluated"
    x::Vector{SVector{Dim,T}}
    "Field values at corresponding positions"
    field::Vector{SVector{FieldDim,T}}

    function FieldResults(x::AbstractVector{<:AbstractVector}, field::AbstractVector{<:AbstractVector})
        T = eltype(x[1])
        Dim = length(x[1])
        FieldDim = length(field[1])

        # Input validation
        if length(field) != length(x)
            throw(ArgumentError("The field and spatial positions must have the same number of elements"))
        end
        if !all(p -> length(p) == Dim, x)
            throw(ArgumentError("All positions must have dimension $Dim"))
        end
        if !all(f -> length(f) == FieldDim, field)
            throw(ArgumentError("All field values must have dimension $FieldDim"))
        end

        new{T,Dim,FieldDim}(convert(Vector{SVector{Dim,T}}, x),
                           convert(Vector{SVector{FieldDim,T}}, field))
    end
end

field(fr::FieldResults) = fr.field

struct FundamentalSolution{T<:Real,Dim,FieldDim}
    source_positions::Vector{SVector{Dim,T}}
    coefficients:: Vector{SVector{FieldDim,T}}
end

# Interface implementations
Base.length(fr::FieldResults) = length(fr.x)
Base.getindex(fr::FieldResults, i::Int) = (fr.x[i], fr.field[i])
Base.iterate(fr::FieldResults) = length(fr) == 0 ? nothing : ((fr.x[1], fr.field[1]), 1)
Base.iterate(fr::FieldResults, state) = state >= length(fr) ? nothing : ((fr.x[state+1], fr.field[state+1]), state + 1)