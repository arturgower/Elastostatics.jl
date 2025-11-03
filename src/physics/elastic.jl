"""
    Elastic{Dim,T<:AbstractFloat}(ρ::T, c::Complex{T})
    Elastic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic elastic medium for static simulations.

Simulations in this medium produce scalar (Dim) fields in Dim dimensions.
"""
struct Elasticstatic{Dim,T} <: PhysicalMedium{Dim,Dim}
    ρ::T # Density (use \greekletter+tab to get the greek letter!)
    cp::Complex{T} # Phase velocity of pressure wave
    cs::Complex{T} # Phase velocity of shear wave
end

# Constructor which supplies the dimension without explicitly mentioning type
function Elasticstatic(Dim::Integer; ρ::T = 0.0, cp::Union{T,Complex{T}} = 0.0, cs::Union{T,Complex{T}} = 0.0) where {T<:Number}
     Elastic{Dim,T}(ρ,Complex{T}(cp),Complex{T}(cs))
end

struct DisplacementType <: FieldType end
struct TractionType <: FieldType end
