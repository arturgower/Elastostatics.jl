"""
    Elastic{Dim,T<:AbstractFloat}(ρ::T, c::Complex{T})
    Elastic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic elastic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (Dim) fields in Dim dimensions. In general we use the Debye potentials to describe the field.
"""
struct Elastic{Dim,T} <: PhysicalMedium{Dim,Dim}
    ρ::T # Density (use \greekletter+tab to get the greek letter!)
    cp::Complex{T} # Phase velocity of pressure wave
    cs::Complex{T} # Phase velocity of shear wave
end

# Constructor which supplies the dimension without explicitly mentioning type
function Elastic(Dim::Integer; ρ::T = 0.0, cp::Union{T,Complex{T}} = 0.0, cs::Union{T,Complex{T}} = 0.0) where {T<:Number}
     Elastic{Dim,T}(ρ,Complex{T}(cp),Complex{T}(cs))
end


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

# Constructor without type parameters
function BoundaryData(field_type::F; 
        points = [zeros(3)], 
        interior_points = mean(points) |> x -> [SVector(x...)],
        fields = [p .* 0.0 for p in points]
    ) where F <: FieldType

    Dim = length(points[1])

    return BoundaryData{F,Dim}(interior_points, points, fields)
end