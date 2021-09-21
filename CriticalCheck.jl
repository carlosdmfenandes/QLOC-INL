include("Permanent.jl")

"Checks if a matrix is a stationary point of the permanent."

"Calculates the permanent adjoint"
function padjoint(mat)
    iter = 1:size(mat,1)
    val = [permanent(mat[1:end .!=i,1:end .!=j])
           for i in iter,
               j in iter
           ]
    return val
end

"Calculates how far mat is form being Hermitian."
function criticalcheck(mat)
    testmat = transpose(mat)*padjoint(mat)
    LinearAlgebra.norm(testmat-adjoint(testmat))
end
