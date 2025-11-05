

function FundamentalSolution(medium::P, bd::BoundaryData{F,Dim}; tol = 1e-5) where {P <: PhysicalMedium, F <: FieldType, Dim}

    using BlockArrays
    # Matrix(mortar(Matrix{Matrix{Complex{T}}}(E1s)))
    positions = source_positions(cloud)
    
    Ms = [
        greens(medium, F, bd.points[i] - x, bd.outward_normals[i])    
    for i in eachindex(bd.points), x in positions]

    forcing = vcat(bd.fields...)

    Matrix(mortar(Ms))

    # δ = method.regularisation_parameter
            # bigA = [A; sqrt(δ) * I];
            # x = bigA \ [b; zeros(size(A)[2])]

            # condition matrix
    SM = diagm([4.0 / sum(abs.(BB[:,j])) for j in 1:size(BB,2)])
    BBSM = BB * SM

    # coes = M \ forcing;

    # Tikinov
    δ = sqrt(tol * norm(forcing) /  maximum(size(A)))
    bigM = [A; sqrt(δ) * I];
    x = bigA \ [forcing; zeros(size(A)[2])]

    mortar(Ms)

end
