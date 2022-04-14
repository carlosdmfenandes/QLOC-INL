struct BeamSplitter <: AbstractMatrix{ComplexF64}
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
size(bs::BeamSplitter) = ntuple(x -> max(bs.m, bs.n), Val(2))

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
#=
"Extend rand to BeamSplitters."
function rand(BeamSplitter, (n, m), transpose=true)
    BeamSplitter(m, n, 2*pi*rand(), 2*pi*rand(), transpose)
end

"Generate random BS acting on random modes. Code is straightfoward
but could be made to call the RNG only once."
function rand(BeamSplitter, M::Integer, transpose=true)
    n = rand(1:M)
    d = rand(1:(M-1))
    m = (n+d)%M
    rand(bs, (n, m), transpose)
end
=#
"Extend matrix multiplication to include Beam Splitters."
*(A::BeamSplitter, B::AbstractMatrix) = x!(A, copy(B))
*(A::AbstractMatrix, B::BeamSplitter) = x!(copy(A), B)

