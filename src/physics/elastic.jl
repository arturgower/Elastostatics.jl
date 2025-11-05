"""
    Elastostatic{Dim,T<:AbstractFloat}(ρ::T, c::Complex{T})
    Elastostatic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic elastic medium for static simulations.

Simulations in this medium produce scalar (Dim) fields in Dim dimensions.
"""
struct Elastostatic{Dim,T} <: PhysicalMedium{Dim,Dim}
    ρ::T # Density (use \greekletter+tab to get the greek letter!)
    cp::Complex{T} # Phase velocity of pressure wave
    cs::Complex{T} # Phase velocity of shear wave
    
    # Constructor which supplies the dimension without explicitly mentioning type
    function Elastostatic(Dim::Integer; ρ::T = 0.0, cp::Union{T,Complex{T}} = 0.0, cs::Union{T,Complex{T}} = 0.0) where {T<:Number}
        # check Lame parameters are positive
        μ = ρ * cs^2  
        λ = ρ * cp^2 - 2μ

        if μ < 0 || λ < 0 
            @error "The Lame parameters need to be positive."
        end

        new{Dim,T}(ρ,Complex{T}(cp),Complex{T}(cs))
    end
end

struct DisplacementType <: FieldType end
struct TractionType <: FieldType end

function greens(medium::Elastostatic{2,T}, FT::TractionType, x::SVector{2,T}, outward_normal::SVector{2,T}) where T

    n = outward_normal
    μ = medium.ρ * medium.cs^2  
    λ = medium.ρ * medium.cp^2 - 2μ

    r2 = dot(x,x)
    xn = dot(x,n)

    Σn = [
        μ * (xn * (l == i) + x[i]*n[l] - x[l]*n[i]) + 2(λ+μ) * (x[l]*x[i]*xn / r2)
    for i = 1:2, l = 1:2]
    Σn = Σn ./ (2pi*(λ+2μ)*r2)

    return Σn 
end

function greens(medium::Elastostatic{2,T}, FT::DisplacementType, x::SVector{2,T}, outward_normal::SVector{2,T}) where T

    μ = medium.ρ * medium.cs^2  
    λ = medium.ρ * medium.cp^2 - 2μ

    ν = λ / (2λ + 2μ) 
    r2 = dot(x,x)

    U = [
        (3-4ν) * log(r) * (i==j) - x[i] * x[j] / r2
    for i = 1:2, j = 1:2]

    U = U ./ (8pi*μ*(1-ν))
    return U
end