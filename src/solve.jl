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
struct FieldResults{Dim,FieldDim,T<:Real}
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

        new{Dim,FieldDim,T}(convert(Vector{SVector{Dim,T}}, x),
                           convert(Vector{SVector{FieldDim,T}}, field))
    end
end

field(fr::FieldResults) = fr.field

# Interface implementations
Base.length(fr::FieldResults) = length(fr.x)
Base.getindex(fr::FieldResults, i::Int) = (fr.x[i], fr.field[i])
Base.iterate(fr::FieldResults) = length(fr) == 0 ? nothing : ((fr.x[1], fr.field[1]), 1)
Base.iterate(fr::FieldResults, state) = state >= length(fr) ? nothing : ((fr.x[state+1], fr.field[state+1]), state + 1)

struct FundamentalSolution{Dim,P<:PhysicalMedium{Dim},T<:Real}
    medium::P
    positions::Vector{SVector{Dim,T}}
    coefficients::Array{T}

    function FundamentalSolution(medium::P,
            positions::Vector{<:AbstractVector},
            coefficients::AbstractArray) where {P<:PhysicalMedium}
        
        # Extract dimension information
        Dim = dimension(medium)
        T = eltype(positions[1])
        
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

        return new{Dim,P,T}(medium, pos_converted, coefficients)
    end
end

function FundamentalSolution(medium::P, bd::BoundaryData{F,Dim}; tol = 1e-5) where {P <: PhysicalMedium, F <: FieldType, Dim}

    positions = source_positions(bd)
    
    Ms = [
        greens(medium, F, bd.boundary_points[i] - x, bd.outward_normals[i])    
    for i in eachindex(bd.boundary_points), x in positions]

    forcing = vcat(bd.fields...)

    Matrix(mortar(Ms))

    # δ = method.regularisation_parameter
            # bigA = [A; sqrt(δ) * I];
            # x = bigA \ [b; zeros(size(A)[2])]

            # condition matrix
    SM = diagm([4.0 / sum(abs.(BB[:,j])) for j in 1:size(BB,2)])
    BBSM = BB * SM

    # coes = M \ forcing;

    # Tikinov
    δ = sqrt(tol * norm(forcing) /  maximum(size(M)))
    bigM = [M; sqrt(δ) * I];
    x = bigA \ [forcing; zeros(size(A)[2])]

    mortar(Ms)

end
