

function FundamentalSolution(bd::BoundaryData{F,Dim}) where {F <: FieldType, Dim}

using BlockArrays
# Matrix(mortar(Matrix{Matrix{Complex{T}}}(E1s)))
Matrix(mortar(Ms))

mortar(Ms)

end
