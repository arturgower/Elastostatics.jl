function interior_points_along_coordinate(points; coordinate_number = 2, segments = 10)
    
    cs = [xy[coordinate_number] for xy in points]

    minc = minimum(cs)
    maxc = maximum(cs)
    edges = collect(range(minc, stop=maxc, length=segments+1))

    # collect indices for each y-bin
    bin_indices = [Int[] for _ in 1:segments]
    
    for (idx, y) in enumerate(cs)
        b = clamp(searchsortedlast(edges, y), 1, segments)
        push!(bin_indices[b], idx)
    end

    interior_points = map(bin_indices) do inds
        mean(points[inds])
    end

    return interior_points
end
