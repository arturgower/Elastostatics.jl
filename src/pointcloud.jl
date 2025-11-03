"""
    PointCloud(boundary_points::Vector{SVector{Dim,T}}, interior_points::Vector{SVector{Dim,T}}) where {T, Dim}

A [`Shape`](@ref) defined by a set of points on the boundary `boundary_points`. There is no particular order to the `boundary_points`. The `interior_points` are a set of points that are in the interior of the body and are used to quickly determine what is inside or not. See
@doc in(x::AbstractVector, cloud::PointCloud).
"""
struct PointCloud{T,Dim} <: Shape{Dim}
    boundary_points::Vector{SVector{Dim,T}}
    outward_normals::Vector{SVector{Dim,T}}
    interior_points::Vector{SVector{Dim,T}}
end

function PointCloud(boundary_points::Vector{V}; 
        interior_points::Vector{V} = [mean(boundary_points)],
        outward_normals::Vector{V} = [boundary_points[1] .* zero(typeof(boundary_points[1][1]))]
    ) where {V <: AbstractVector}

    T = eltype(boundary_points[1])
    Dim = length(boundary_points[1])
    Polygon{T,Dim}(boundary_points, outward_normals, interior_points)
end

import MultipleScattering: name
name(shape::PointCloud) = "PointCloud"

bounding_box(poly::PointCloud) = Box(poly.boundary_points)

# import Base.in
# function in(x::AbstractVector, poly::Polygon)::Bool
#     # Ray-casting algorithm for point-in-polygon test
#     n = length(poly.boundary_points)
#     inside = false
#     j = n
#     for i in 1:n
#         xi, yi = poly.boundary_points[i][1], poly.boundary_points[i][2]
#         xj, yj = poly.boundary_points[j][1], poly.boundary_points[j][2]
#         if ((yi > x[2]) != (yj > x[2])) && (x[1] < (xj - xi) * (x[2] - yi) / (yj - yi) + xi)
#             inside = !inside
#         end
#         j = i
#     end
#     return inside
# end

import Base.issubset
function issubset(poly::PointCloud, box::Box)
    poly_box = bounding_box(poly)
    return issubset(poly_box,box)
end

"""
    issubset(box::Box, poly::Polygon)

Returns true if the corners of the box are contained within polygon, false otherwise.
"""
function issubset(box::Box, poly::PointCloud)
    return all(c ∈ poly for c in corners(box))
end


function outward_normals(cloud::PointCloud)

    pts = cloud.boundary_points
    n = length(pts)
    if n == 0
        return Vector{typeof(pts[1])}()
    end
    Dim = length(pts[1])
    T = eltype(pts[1])

    # number of neighbors to use (at least Dim, at most n-1)
    k = min(Dim+2, max(1, n-1))

    normals = Vector{typeof(pts[1])}(undef, n)

    # helper to compute squared distances
    sqrddist(a,b) = sum((a .- b).^2)

    # for orientation choose the closest interior point to each boundary point
    for i in 1:n
        p = pts[i]

        # find k nearest neighbours (including the point itself)
        dists = [sqrddist(p, q) for q in pts]
        idx = sortperm(dists)[1:min(k+1, n)]   # +1 because p itself is at distance 0
        neighbors = pts[idx]

        # center data and compute covariance
        m = length(neighbors)
        μ = reduce(+, neighbors) / m
        X = hcat((neighbors .- μ)...)            # Dim x m matrix
        C = (X * transpose(X)) / max(1, m-1)    # covariance-like matrix (Dim x Dim)

        # eigen-decomposition: smallest eigenvalue eigenvector is normal to a (Dim-1)-manifold
        E = eigen(Symmetric(C))
        jmin = argmin(E.values)
        nvec = E.vectors[:, jmin]
        if norm(nvec) == 0
            # fallback: use vector from point to centroid of boundary points
            centroid = reduce(+, pts) / n
            nvec = p .- centroid
            if norm(nvec) == 0
                # fallback 2: unit x-axis
                nvec = zeros(Dim); nvec[1] = 1.0
            end
        end
        nvec = nvec / norm(nvec)

        # orient normal to point outward (away from interior)
        # pick closest interior point
        intdists = [sqrddist(p, q) for q in cloud.interior_points]
        intp = cloud.interior_points[argmin(intdists)]
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