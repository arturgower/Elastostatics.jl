using MethodOfFundamentalSolutions

    θs = θs_arr[i]
    r = 1.3
    r = 1.0
    points = [[r*cos(θ), r*sin(θ)] for θ in θs]
    normals = [[cos(θ), sin(θ)] for θ in θs]

    predict_fields = [field(TractionType(), fsol, points[i], normals[i]) for i in eachindex(points)]
    fields = [radial_to_cartesian_transform([r,θ])*[σrr(r,θ), σrθ(r,θ)] for θ in θs]

    errors = [norm(fields[i] - predict_fields[i]) for i in eachindex(fields)]
    errors = errors ./ mean([norm(f) for f in fields])
    
    maximum(errors)

    using Plots
    x_vec, inds = points_in_shape(bd; res = 15)
    x_vec, inds = points_in_shape(bd; res = 25)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.
    xs = x_vec[inds]

    fs = [
        field(TractionType(), fsol, x, x / norm(x)) 
    for x in xs];
    
    field_mat = [[0.0, 0.0] for x in x_vec]
    field_mat[inds] = fs;
    field_predict = FieldResult(x_vec, field_mat[:]);

     # fs = [field(wave,x,fieldtype) for x in xs];
    fs = map(xs) do x
        r, θ = cartesian_to_radial_coordinates(x)
       radial_to_cartesian_transform([r,θ])*[σrr(r,θ), σrθ(r,θ)]
    end

    field_mat = [[0.0, 0.0] for x in x_vec]
    field_mat[inds] = fs;
    field_true = FieldResult(x_vec, field_mat[:]);

    using Plots
    plot(field_predict, field_apply = first, c = :inferno)
    plot(field_true, field_apply = first, c = :inferno)
    
    plot(field_predict - field_true, field_apply = norm, c = :inferno)
    plot!(fsol, xlims = (-1.6,1.6), ylims = (-1.6,1.6))
    plot!(bd)
    
    M = system_matrix(fsol, bd)
    cond(M)
    svdM = svd(M)
    svdM.S[end-6:end]

    norm(M * svdM.Vt[end,:] - svdM.S[end] .* svdM.U[:,end]) / norm(M * svdM.Vt[end,:])

    fsol_null = FundamentalSolution(fsol.medium, fsol.positions, svdM.Vt[end,:])
    
    fs = [
        field(TractionType(), fsol_null, x, x / norm(x)) 
    for x in xs];
    
    field_mat = [[0.0, 0.0] for x in x_vec]
    field_mat[inds] = fs;
    field_null = FieldResult(x_vec, field_mat[:]);

    plot(field_null, field_apply = norm, c = :inferno)