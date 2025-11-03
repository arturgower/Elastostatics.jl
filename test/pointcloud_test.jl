
@testset "point cloud" begin

    # Use points sampled on the boundary with different resolutions
    n = 15
    θs = LinRange(0,2pi,n)[1:(n-1)];

    r = 1.3
    ε = 0.2

    # Create point cloud with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_ns = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        cloud = PointCloud(points; interior_points = interior_points)

        errors = [norm(outward_ns[i] - cloud.outward_normals[i]) for i in eachindex(outward_normals)]
        maximum(errors)

    ε = 0.2
    # Create point cloud with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_ns = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        cloud = PointCloud(points; interior_points = interior_points)

        errors = [norm(outward_ns[i] - cloud.outward_normals[i]) for i in eachindex(outward_normals)]
        maximum(errors)
    
end