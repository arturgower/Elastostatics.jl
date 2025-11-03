
@testset "point cloud" begin

    # Use points sampled on the boundary with different resolutions
    n = 15
    θs = LinRange(0,2pi,n)[1:(n-1)];

    r = 1.3
    ε = 0.2
    ε = 0.0

    # Create point cloud with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_normals = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        cloud = PointCloud(points; outward_normals = [[0.0,0.0]], interior_points = interior_points)

    
end