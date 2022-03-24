module QLOC

using LinearAlgebra

import Base: *, transpose, inv, conj

export BeamSplitter

struct BeamSplitter <: AbstractSparseMatrix{ComplexF64, Int}
    m::Int
    n::Int
    phase::Float64
    rangle::Float64
    transposed::Bool
end

BeamSplitter(m, n, phase, rangle) = BeamSplitter(m, n, phase, rangle, false)
BeamSplitter(m, n) = BeamSplitter(m, n, 0, 0, false)
BeamSplitter() = BeamSplitter(1, 2, 0, 0, false)

"Define standard operations on a beam splitter."
transpose(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, bs.phase, bs.rangle, !bs.transposed)
inv(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, -bs.phase, bs.rangle, !bs.transposed)
conj(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, -bs.phase, bs.rangle, bs.transposed)


"Multiply a matrix by a beamsplitter to the left in place."
function x!(mat::AbstractVecOrMat{ComplexF64}, bs::BeamSplitter)
    ph = cis(bs.phase)
    s, c = sincos(bs.rangle)
    X = view(mat, :, bs.m)
    Y = view(mat, :, bs.n)
    oldY = copy(Y)
    if bs.transposed
        LinearAlgebra.axpby!(s * ph, X, c, Y)
        LinearAlgebra.axpby!(-s, oldY, c * ph, X)
    else
        LinearAlgebra.axpby!(-s, X, c, Y)
        LinearAlgebra.axpby!(s * ph, oldY, c * ph, X)
    end
    return mat
end

"Multiply a matrix by a beamsplitter to the right in place."
function x!(bs::BeamSplitter, mat::AbstractVecOrMat)
    cb = transpose(bs)
    x!(Transpose(mat), cb)
    return mat
end

"Extend matrix multiplication to include Beam Splitters."
*(A::AbstractArray, B::AbstractArray) = A * B
*(A::BeamSplitter, B::AbstractArray) = x!(A, deepcopy(B))
*(A::AbstractArray, B::BeamSplitter) = x!(deepcopy(A), B)

end
