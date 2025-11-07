@testset "boundarydata" begin

    # Use points sampled on the boundary with different resolutions
    n = 15
    θs = LinRange(0,2pi,n)[1:(n-1)];

    r = 1.3
    ε = 0.0

    # Create point cloud with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_ns = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]

        outward_ns2 = outward_normals(points, interior_points)
        errors = [norm(outward_ns[i] - outward_ns2[i]) for i in eachindex(outward_ns)]

    @test maximum(errors) < 1e-10

    ε = 0.1
    # Create point cloud with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_ns = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        
        outward_ns2 = outward_normals(points, interior_points)
        errors = [norm(outward_ns[i] - outward_ns2[i]) for i in eachindex(outward_ns)]

    @test maximum(errors) < 0.1
   
    # Let us check whether the interior is correctly defined.
    cloud = BoundaryData(DisplacementType(); 
        boundary_points = points,
        interior_points = interior_points
    )

    x_vec, inds = points_in_shape(cloud; res = 20)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.

    xs = x_vec[inds]

    @test all([norm(x) < r * 1.1 for x in xs])

    field_mat = [[0.0, 0.0] for x in x_vec]

    # fs = [field(wave,x,fieldtype) for x in xs];
    fs = [[1.0, 1.0] for x in xs];
    field_mat[inds] = fs

    res = FieldResult(x_vec, field_mat[:])

    # using Plots
    # plot(res, clims = (-1,1))
    # plot!(cloud)
end