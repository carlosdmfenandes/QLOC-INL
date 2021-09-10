include("BetterPermanent.jl")

function padjoint(mat)
    iter = 1:size(mat,1)
    val = [permanent(mat[1:end .!=i,1:end .!=j])
           for i in iter,
               j in iter
           ]
    return val
end

function vpadjoint(mat)
    iter = 1:size(mat,1)
    val = [permanent(mat[1:end .!=i,1:end .!=j])
           for i in iter,
               j in iter
           ]
    return val
end

function criticalcheck(mat)
    testmat = transpose(mat)*padjoint(mat)
    LinearAlgebra.norm(testmat-adjoint(testmat))
1end

function criticalcheck(mat::AbstractArray{<:Integer})
    testmat = transpose(mat)*padjoint(mat)
    ishermitian(testmat)
end
