# Here we use solutions from Airy stress function to generate boundary data, and then solve with this package to compare.
@testset "Benchmark circular domain" begin

    medium = Elastostatic(2; ρ = 1.0, cp = 2.0, cs = 1.0)

    # Use points sampled on the boundary with different resolutions
    θs_arr = [LinRange(0,2pi,n)[1:(n-1)] for n in (10,20,40,80)];

    # Circumferential stress
    # From the Airy stress function we have the solution inside a circular domain which does not depend on the radius r: 
    σrr(r, θ) = -2 * cos(2θ)
    σθθ(r, θ) = 2 * cos(2θ)
    σrθ(r, θ) = 2 * sin(2θ)

    r = 1.3

    # Create boundary data for each resolution
    bds = map(θs_arr) do θs
        points = [[r*cos(θ), r*sin(θ)] for θ in θs]
        outward_normals = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]

        # WRRRROOONG needs the basis vectors
        fields = [radial_to_cartesian_transform(SVector(r,θ))*[σrr(r,θ), σrθ(r,θ)] for θ in θs]
        BoundaryData(TractionType(); 
            boundary_points = points, 
            fields = fields, 
            outward_normals = outward_normals,
            interior_points = interior_points
        )
    end

    # Let's check one of the solutions and do lots of plots
    i = 4;
    bd = bds[i]
    fsol = FundamentalSolution(medium, bd; tol = 1e-6);

    using Plots

    x_vec, inds = points_in_shape(bd; res = 15)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.
    xs = x_vec[inds]

    fs = [
        field(TractionType(), fsol, x, x / norm(x)) 
    for x in xs];
    
    field_mat = [[0.0, 0.0] for x in x_vec]
    field_mat[inds] = fs;
    field_predict = FieldResults(x_vec, field_mat[:]);

     # fs = [field(wave,x,fieldtype) for x in xs];
    fs = map(xs) do x
        r, θ = cartesian_to_radial_coordinates(SVector(x...))
       radial_to_cartesian_transform(SVector(r,θ))*[σrr(r,θ), σrθ(r,θ)]
    end

    field_mat = [[0.0, 0.0] for x in x_vec]
    field_mat[inds] = fs;
    field_true = FieldResults(x_vec, field_mat[:]);

    # using Plots
    plot(field_predict, clims = (-1,1), field_apply = first)
    
    plot(field_true, clims = (-1,1), field_apply = first)
    plot!(bd)
    plot!(fsol, xlims = (-1.6,1.6), ylims = (-1.6,1.6))

    θs = θs_arr[i]
    r = 1.3
    r = 1.0
    points = [[r*cos(θ), r*sin(θ)] for θ in θs]
    normals = [[cos(θ), sin(θ)] for θ in θs]

    predict_fields = [field(TractionType(), fsol, points[i], normals[i]) for i in eachindex(points)]
    fields = [radial_to_cartesian_transform(SVector(r,θ))*[σrr(r,θ), σrθ(r,θ)] for θ in θs]

    errors = [norm(fields[i] - predict_fields[i]) / norm(fields[i]) for i in eachindex(fields)]
    maximum(errors)
    

    # From the Airy stress function we another the solution inside a circular domain which does depend on the radius r: 
    σrr(r,θ) = 0.0;
    σθθ(r,θ) = 12 * r^2 * cos(2θ)
    σrθ(r,θ) = 6 * r^2 * sin(2θ)

    r = 1.3

    # Create boundary data for each resolution
    bds = map(θs_arr) do θs
        points = [[r*cos(θ), r*sin(θ)] for θ in θs]
        outward_normals = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        fields = [radial_to_cartesian_transform(SVector(r,θ))*[σrr(r,θ), σrθ(r,θ)] for θ in θs]
        BoundaryData(TractionType(); 
            boundary_points = points, 
            fields = fields, 
            outward_normals = outward_normals,
            interior_points = interior_points
        )
    end

    i = 4;
    bd = bds[i]
    fsol = FundamentalSolution(medium, bd; tol = 1e-6)
    
    θs = θs_arr[i]
    r = 1.3
    r = 0.3
    points = [[r*cos(θ), r*sin(θ)] for θ in θs]
    normals = [[cos(θ), sin(θ)] for θ in θs]

    predict_fields = [field(TractionType(), fsol, points[i], normals[i]) for i in eachindex(points)]
    fields = [radial_to_cartesian_transform(SVector(r,θ))*[σrr(r,θ), σrθ(r,θ)] for θ in θs]

    errors = [norm(fields[i] - predict_fields[i]) / norm(fields[i]) for i in eachindex(fields)]
    maximum(errors)
    

    
end

@testset "Benchmark annulus" begin
    # From the Airy stress function we have the solution inside an annulus domain: 
    σrr(r,θ) = -2*cos(θ) / r^3
    σθθ(r,θ) = 2*cos(θ) / r^3
    σrθ(r,θ) = - 2*sin(θ) / r^3
end

