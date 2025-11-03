@testset "point cloud" begin

    # Use points sampled on the boundary with different resolutions
    n = 15
    θs = LinRange(0,2pi,n)[1:(n-1)];

    r = 1.3
    ε = 0.0

    # Create point cloud with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_ns = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        cloud = PointCloud(points; interior_points = interior_points)

        errors = [norm(outward_ns[i] - cloud.outward_normals[i]) for i in eachindex(outward_ns)]
    @test maximum(errors) < 1e-10

    ε = 0.1
    # Create point cloud with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_ns = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        cloud = PointCloud(points; interior_points = interior_points)

        #using Plots
        # plot(cloud)

        errors = [norm(outward_ns[i] - cloud.outward_normals[i]) for i in eachindex(outward_ns)]

    @test maximum(errors) < 0.1
   
    
    x_vec, inds = points_in_shape(cloud; res = 12)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.

    xs = x_vec[inds]
    field_mat = [SVector(0.0+0.0im, 0.0+0.0im) for x in x_vec]

    fs = [field(wave,x,fieldtype) for x in xs];
    field_mat[inds] = fs

    return  FrequencySimulationResult(reshape(field_mat, :, 1), x_vec, [wave.ω])
end