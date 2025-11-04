# recipe for plotting a PointCloud. Assumes 2D boundary_points 
@recipe function plot(cloud::PointCloud)
    # Set default attributes
    markershape --> :circle
    legend --> false
    aspect_ratio --> 1.0

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
        markersize --> 4.0
        (bx, by)
    end

    @series begin
        label --> "Interior points"
        seriestype --> :scatter
        markersize --> 6.0
        (ix, iy)
    end

    @series begin
        label --> ""
        color --> :black
        seriestype --> :quiver
        quiver --> (nx, ny)
        marker = :none           # remove marker symbol at arrow tips
        markersize --> 0.0          # fallback to ensure no marker is drawn
        (bx, by)
    end
end


# Plot the result in space (across all x) for a specific angular frequency
@recipe function plot(res::FieldResults;
        region_shape = :empty, field_apply := first)

    x = [x[1] for x in res.x]
    y = [x[end] for x in res.x] # y will actually be z for 3D...

    seriestype --> :heatmap
    seriescolor --> :balance
    aspect_ratio --> 1.0

    st = get(plotattributes, :seriestype, :surface)
    
    if st == :heatmap
        # We could check here to see if x and y have the right structure
        x = unique(x)
        y = unique(y)

        n_x = length(x)
        n_y = length(y)

        fill --> true

        if region_shape != :empty
            bounds = bounding_box(region_shape)

            # If user has not set xlims and ylims, set them to the rectangle
            xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
            ylims --> (bottomleft(bounds)[2], topright(bounds)[2])
        else
            xlims --> (minimum(x), maximum(x))
            ylims --> (minimum(y), maximum(y))
        end

        x, y, field_apply.(transpose(reshape(field(res),n_x,n_y)))

    else
        (x, y, field_apply.(field(res)))
    end

end