include("../HadamardMatrices.jl")
include("../QLOC.jl")

using .QLOC

const tol = 1e-15
const sillym = Complex.(sylvester_matrix(2)/2)
const ox = OxfordDecomp(sillym)

@show abs2.(Matrix(ox) - sillym) |> sum <= tol
qtikzbody(ox.bsvec)
