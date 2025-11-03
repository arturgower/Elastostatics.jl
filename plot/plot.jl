using RecipesBase

# recipe for plotting a PointCloud. Assumes 2D boundary_points 
@recipe function plot(cloud::PointCloud)

    bs = cloud.boundary_points
    ns = cloud.outward_normals
    is = cloud.interior_points

    # used to determine the length to plot the vectors
    lengthscale = mean([norm(b - mean(bs)) for b in bs])

    bx = [b[1] for b in bs]
    by = [b[2] for b in bs]
    ix = [b[1] for b in is]
    iy = [b[2] for b in is]

    nx = [b[1] for b in ns]
    ny = [b[2] for b in ns]

    nx = nx .* (lengthscale/4)
    ny = ny .* (lengthscale/4)

    # fallback: repeat or zero-length normals
    if length(nx) != length(bx)
        nx = zeros(length(bx))
        ny = zeros(length(by))
    end

    @series begin
        label --> "Boundary points"
        seriestype --> :scatter
        markersize -> 2
        (bx, by)
    end

    @series begin
        label --> ""
        color --> :black
        seriestype --> :quiver
        quiver --> (nx, ny)
        (bx, by)
    end

    @series begin
        label --> "Interior points"
        seriestype --> :scatter
        markersize -> 4
        (ix, iy)
    end
end

"Plot the field for a particular wavenumber"
@recipe function plot(res::Simulation{2};
        resolution = 10, res = resolution, xres=res, yres=res,
        field_apply=real,
        region_shape = :auto,
        bounds = :auto,
        exclude_region = EmptyShape{2}(),
        drawparticles=true)

    # If user wants us to, generate bounding rectangle around particles
    region_shape = (region_shape != :auto) ? region_shape :
        if isempty(sim.particles)
            if bounds == :auto
                @warn "What region to plot? For example, use keyword bounds = Box([[-1.0,-1.0],[1.0,1.0]])"
                Box([[-1.0,-1.0],[1.0,1.0]])
            else bounds
            end
        else
            region_shape = bounding_box(sim.particles)
        end

    bounds = bounding_box(region_shape)
    # If user has not set xlims and ylims, set them to the rectangle
    xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
    ylims --> (bottomleft(bounds)[2], topright(bounds)[2])

    # Incase the user did set the xlims and ylims, generate a new bounding
    # rectangle with them
    p_xlims = plotattributes[:xlims]
    p_ylims = plotattributes[:ylims]
    bounds = Box([[p_xlims[1],p_ylims[1]], [p_xlims[2],p_ylims[2]]])

    region_shape = (bounds ⊆ region_shape) ? bounds : region_shape

    field_sim = run(sim, region_shape, [ω]; xres=xres, yres=yres, zres=yres, exclude_region=exclude_region)
    xy_mat = reshape(field_sim.x, (xres+1, yres+1))
    x_pixels = [x[1] for x in xy_mat[:,1]]
    y_pixels = [x[2] for x in xy_mat[1,:]]

    @series begin

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field(field_sim), (xres+1, yres+1)))
        seriestype --> :contour
        fill --> true
        grid --> false
        aspect_ratio := 1.0
        seriescolor --> :balance
        title --> "Field at ω=$ω"

        (x_pixels, y_pixels, field_apply.(response_mat))
    end

    if drawparticles
        @series begin
            sim.particles
        end
    end

end