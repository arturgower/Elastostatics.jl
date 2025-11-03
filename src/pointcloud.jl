"""
    PointCloud(boundary_points::Vector{SVector{Dim,T}}, interior_points::Vector{SVector{Dim,T}}) where {T, Dim}

A [`Shape`](@ref) defined by a set of points on the boundary `boundary_points`. There is no particular order to the `boundary_points`. The `interior_points` are a set of points that are in the interior of the body and are used to quickly determine what is inside or not. See
@doc in(x::AbstractVector, cloud::PointCloud).
"""
struct PointCloud{T,Dim} <: Shape{Dim}
    boundary_points::Vector{SVector{Dim,T}}
    outward_normals::Vector{SVector{Dim,T}}
    interior_points::Vector{SVector{Dim,T}}

    function PointCloud(boundary_points::Vector{V}; 
            interior_points::Vector{V} = [mean(boundary_points)],
            outward_normals::Vector{V} = outward_normals(boundary_points, interior_points)
        ) where {V <: AbstractVector}

        T = eltype(boundary_points[1])
        Dim = length(boundary_points[1])

        new{T,Dim}(boundary_points, outward_normals, interior_points)
    end
end

# function PointCloud(boundary_points::Vector{V}; 
#         interior_points::Vector{V} = [mean(boundary_points)],
#         outward_normals::Vector{V} = [boundary_points[1] .* zero(typeof(boundary_points[1][1]))]
#     ) where {V <: AbstractVector}

#     T = eltype(boundary_points[1])
#     Dim = length(boundary_points[1])
#     PointCloud{T,Dim}(boundary_points, outward_normals, interior_points)
# end


import MultipleScattering: name
name(shape::PointCloud) = "PointCloud"

bounding_box(cloud::PointCloud) = Box(cloud.boundary_points)

import Base.in
function in(x::AbstractVector, cloud::PointCloud)::Bool

    # find nearest point p on the boundary to x. And nearest interior point q to x. If |q - x| < |q - p| then the point is inside the body defined by cloud 
    dists = [sum((x - p) .^2) for p in cloud.boundary_points];
    p = cloud.boundary_points[argmin(dists)]

    dists = [sum((x - p) .^2) for p in cloud.interior_points];
    q = cloud.interior_points[argmin(dists)]
    
    inside = norm(x - q) < norm(p - q) ? true : false
    return inside
end

import Base.issubset
function issubset(cloud::PointCloud, box::Box)
    cloud_box = bounding_box(cloud)
    return issubset(cloud_box,box)
end

"""
    issubset(box::Box, poly::PointCloud)

Returns true if the corners of the box are contained within polygon, false otherwise.
"""
function issubset(box::Box, cloud::PointCloud)
    return all(c ∈ cloud for c in corners(box))
end


outward_normals(cloud::PointCloud) = outward_normals(cloud.boundary_points, cloud.interior_points)

function outward_normals(boundary_points, interior_points)

    pts = boundary_points
    n = length(pts)
    if n == 0
        return Vector{typeof(pts[1])}()
    end
    Dim = length(pts[1])
    T = eltype(pts[1])

    # number of neighbors to use (at least Dim+1, at most n-1)
    k = min(2Dim, max(1, n-1))

    normals = Vector{typeof(pts[1])}(undef, n)

    # for orientation choose the closest interior point to each boundary point
    for i in 1:n
        p = pts[i]

        # find k nearest neighbours (including the point itself)
        dists = [sum((p - q) .^2) for q in pts]
        idx = sortperm(dists)[1:min(k+1, n)]   # +1 because p itself is at distance 0
        neighbors = pts[idx]

        # center data and compute covariance
        m = length(neighbors)
        μ = mean(neighbors)
        X = hcat([μ - q for q in neighbors]...) # Dim x k+1 matrix
        C = (X * transpose(X)) / max(1, m-1)   # covariance-like matrix (Dim x Dim)

        # eigen-decomposition: smallest eigenvalue eigenvector is normal to a (Dim-1)-manifold
        E = eigen(C)
        jmin = argmin(E.values)
        nvec = E.vectors[:, jmin]
        nvec = nvec / norm(nvec)

        # orient normal to point outward (away from interior)
        # pick closest interior point
        intdists = [sum((p - q) .^2) for q in interior_points]
        intp = interior_points[argmin(intdists)]
        interior_vector = intp .- p    # points from boundary point into interior
        # if dot(nvec, interior_vector) > 0 then nvec points inward -> flip
        if dot(nvec, interior_vector) > 0
            nvec = -nvec
        end

        # convert to same element type as pts
        normals[i] = convert(SVector{Dim,T}, nvec)
    end

    return normals
end