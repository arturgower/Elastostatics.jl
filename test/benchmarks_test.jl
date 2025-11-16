@testset "Benchmark analytic solutions" begin

# Here we use solutions from Airy stress function to generate boundary data, and then solve with this package to compare.
@testset "circular domain" begin

    medium = Elastostatic(2; ρ = 1.0, cp = 2.0, cs = 1.0)

    # Use points sampled on the boundary with different resolutions
    θs_arr = [LinRange(0,2pi,n)[1:(n-1)] for n in (10,80)];

    # Circumferential stress
    # From the Airy stress function we have the solution inside a circular domain which does not depend on the radius r: 
    σrr(r, θ) = -2 * cos(2θ)
    σθθ(r, θ) = 2 * cos(2θ)
    σrθ(r, θ) = 2 * sin(2θ)

    r = 1.3

    # Create boundary data for each resolution
    bds = map(θs_arr) do θs
        
        points = [[r*cos(θ), r*sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]

        normals = map(θs) do θ 
            outward_normal = [cos(θ), sin(θ)] 
        end
        
        fields = map(θs) do θ 
            radial_to_cartesian_transform([r,θ])*[σrr(r,θ), σrθ(r,θ)]
        end

        BoundaryData(TractionType();
            boundary_points = points,
            fields = fields, 
            outward_normals = normals, 
            interior_points = interior_points
        )
    end

    # Solve
    solver = TikhonovSolver(tolerance = 1e-12) 
    
    problems = map(bds) do bd
        Simulation(medium, bd; 
            solver = solver,
            source_positions = source_positions(bd; relative_source_distance = 2.0) 
        )
    end
    fsols = solve.(problems)

    # Predict 
    errors = map(eachindex(fsols)) do n
        θs = θs_arr[n]
        r = 1.0
        points = [[r*cos(θ), r*sin(θ)] for θ in θs]
        normals = [[cos(θ), sin(θ)] for θ in θs]

        predict_fields = [
            field(TractionType(), fsols[n], points[i], normals[i]) 
        for i in eachindex(points)]

        fields = [radial_to_cartesian_transform([r,θ])*[σrr(r,θ), σrθ(r,θ)] for θ in θs]
        error = [norm(fields[i] - predict_fields[i]) for i in eachindex(fields)]
        error ./ mean([norm(f) for f in fields])
        maximum(error)
    end
    
    conditions = map(problems) do p
        system_matrix(p) |> cond
    end
    
    @test errors[1] < 0.1
    @test errors[end] < 0.001

    # From the Airy stress function we have another the solution inside a circular domain which does depend on the radius r: 
    σrr(r,θ) = 0.0;
    σθθ(r,θ) = 12 * r^2 * cos(2θ)
    σrθ(r,θ) = 6 * r^2 * sin(2θ)

    # test adding a particular solution to the boundary data
    
    bd = bds[1]
    psol = ParticularGravity()
    # bd_with_particular = BoundaryData(medium, bd, psol) 
    @test true
end

@testset "Gravity rectangle" begin
    @test true
end

@testset "Circular annulus" begin
    # From the Airy stress function we have the solution inside an annulus domain: 
    σrr(r,θ) = -2*cos(θ) / r^3
    σθθ(r,θ) = 2*cos(θ) / r^3
    σrθ(r,θ) = - 2*sin(θ) / r^3
end

end