# Here we use solutions from Airy stress function to generate boundary data, and then solve with this package to compare.
@testset "Benchmark circular domain" begin

    # Use points sampled on the boundary with different resolutions
    θs_arr = [LinRange(0,2pi,n)[1:(n-1)] for n in (10,20,40,80)];

    # Circumferential stress
    # From the Airy stress function we have the solution inside a circular domain which does not depend on the radius r: 
    σrr(r, θ) = -2*cos(2θ)
    σθθ(r, θ) = 2 * cos(2θ)
    σrθ(r, θ) = 2 * sin(2θ)

    r = 1.3

    # Create point cloud each resolution
    clouds = map(θs_arr) do θs
        points = [[r*cos(θ), r*sin(θ)] for θ in θs]
        outward_normals = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        PointCloud(points; outward_normals = outward_normals, interior_points = interior_points)
    end

    # Create boundary data for each resolution
    bds = map(θs_arr) do θs
        points = [[r*cos(θ), r*sin(θ)] for θ in θs]
        fields = [[σrr(r,θ), σrθ(r,θ)] for θ in θs]
        BoundaryData(TractionType(); points = points, fields = fields)
    end

    bd = bds[1]



    # From the Airy stress function we another the solution inside a circular domain which does depend on the radius r: 
    σrr(r,θ) = 0.0;
    σθθ(r,θ) = 12 * r^2 * cos(2θ)
    σrθ(r,θ) = 6 * r^2 * sin(2θ)

    
end

@testset "Benchmark annulus" begin
    # From the Airy stress function we have the solution inside an annulus domain: 
    σrr(r,θ) = -2*cos(θ) / r^3
    σθθ(r,θ) = 2*cos(θ) / r^3
    σrθ(r,θ) = - 2*sin(θ) / r^3
end

