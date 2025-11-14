
"""
    FundamentalSolution{Dim,P<:PhysicalMedium{Dim},T,C}

A type representing a fundamental solution. The fundamental solution is constructed by placing source points at specified positions outside the domain of interest and determining their coefficients to satisfy boundary conditions. in terms of point sources placed outside the body.

# Fields
- `medium::P`: Physical medium containing material properties
- `positions::Vector{SVector{Dim,T}}`: Vector of positions of the point sources in `Dim`-dimensional space
- `coefficients::Vector{C}`: Vector of coefficients for the point sources
"""
struct FundamentalSolution{Dim,P<:PhysicalMedium{Dim},T,C}
    medium::P
    positions::Vector{SVector{Dim,T}}
    coefficients::Vector{C}

    function FundamentalSolution(medium::P,
            positions::Vector{<:AbstractVector},
            coefficients::AbstractVector) where {P<:PhysicalMedium}
        
        # Extract dimension information
        Dim = spatial_dimension(medium)
        T = eltype(positions[1])
        C = eltype(coefficients)
        
        # Validate dimensions
        if !all(p -> length(p) == Dim, positions)
            throw(ArgumentError("All positions must have dimension $Dim"))
        end

        # Validate coefficients
        FD = field_dimension(medium)
        if length(coefficients) != length(positions) * FD 
            throw(ArgumentError(
                "Expected $(length(positions) * FD) coefficients but got $(length(coefficients))"
            ))
        end

        # Convert inputs to proper types
        pos_converted = convert(Vector{SVector{Dim,T}}, positions)
        coef_converted = convert(Vector{C}, coefficients)

        return new{Dim,P,T,C}(medium, pos_converted, coef_converted)
    end
end

"""
    FieldResult{T<:Real,Dim,FieldDim}

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
results = FieldResult(positions, fields)
```
"""
struct FieldResult{Dim,FieldDim,T<:Real}
    "Spatial positions where field is evaluated"
    x::Vector{SVector{Dim,T}}
    "Field values at corresponding positions"
    field::Vector{SVector{FieldDim,T}}

    function FieldResult(x::AbstractVector{<:AbstractVector}, field::AbstractVector{<:AbstractVector})
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

        new{Dim,FieldDim,T}(convert(Vector{SVector{Dim,T}}, x),
                           convert(Vector{SVector{FieldDim,T}}, field))
    end
end

field(fr::FieldResult) = fr.field

# Interface implementations
import Base.(+)
import Base.(-)
function +(fr1::FieldResult{D,F}, fr2::FieldResult{D,F}) where {D,F}
    # lengths must match
    if length(fr1.x) != length(fr2.x)
        throw(ArgumentError("FieldResult must have the same number of positions"))
    end

    # positions must match elementwise
    if !all(t -> t[1] == t[2], zip(fr1.x, fr2.x))
        throw(ArgumentError("Spatial positions `x` are not compatible between the two FieldResult"))
    end

    # subtract field values elementwise
    newfields = [fr1.field[i] + fr2.field[i] for i in eachindex(fr1.field)]

    return FieldResult(fr1.x, newfields)
end
-(fr1::FieldResult, fr2::FieldResult) = fr1 + FieldResult(fr2.x, - fr2.field)

Base.length(fr::FieldResult) = length(fr.x)
Base.getindex(fr::FieldResult, i::Int) = (fr.x[i], fr.field[i])
Base.iterate(fr::FieldResult) = length(fr) == 0 ? nothing : ((fr.x[1], fr.field[1]), 1)
Base.iterate(fr::FieldResult, state) = state >= length(fr) ? nothing : ((fr.x[state+1], fr.field[state+1]), state + 1)
