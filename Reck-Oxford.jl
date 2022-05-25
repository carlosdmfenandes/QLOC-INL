
@warn """

    DIRECT USE OF "Reck-Oxford.jl" IS DEPRECATED.
    USE QLOC MODULE INSTEAD.
    """

using LinearAlgebra

"""
```
struct BeamSplitter
    m::Int
    n::Int
    phase::Float64
    rangle::Float64
    transposed::Bool

end
```

A type representing a BeamSplitter together with a phase shifter.
The first two fields `m` and `n` represent the modes coupled by the
Beamsplitter.

"""
struct BeamSplitter
    m::Int
    n::Int
    phase::Float64
    rangle::Float64
    transposed::Bool
end

BeamSplitter(m, n, phase, rangle) = BeamSplitter(m, n, phase, rangle, false)
BeamSplitter(m, n) = BeamSplitter(m, n, 0, 0, false)
BeamSplitter() = BeamSplitter(1, 2, 0, 0, false)

"Define standard matrix operations on a beam splitter."
transpose(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, bs.phase, bs.rangle, !bs.transposed)
inverse(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, -bs.phase, bs.rangle, !bs.transposed)
conjugate(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, -bs.phase, bs.rangle, bs.transposed)

"Multiply a matrix by a beamsplitter to the right in place."
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

"Multiply a matrix by a beamsplitter to the left in place."
function x!(bs::BeamSplitter, mat::AbstractVecOrMat)
    cb = transpose(bs)
    x!(Transpose(mat), cb)
    return mat
end

"Extend matrix multiplication to include Beam Splitters."
⊙(A::AbstractArray, B::AbstractArray) = A * B
⊙(A::BeamSplitter, B::AbstractArray) = x!(A, deepcopy(B))
⊙(A::AbstractArray, B::BeamSplitter) = x!(deepcopy(A), B)

"""
Eliminate the matrix element at coordinates 'n, m' by multiplying a
beam splitter.
The beamsplitter pivots on line k if multiplied at the left or
at column k if multiplied at the right.
"""
function elim(mat, n, m, k, left=true)
    dim = size(mat)[1]
    c = mat[n, m]
    if left
        a = mat[k, m]
        rang = atan(abs(c) , abs(a))
        ph = abs2(a) <= 0.0 ? 0.0 : -angle(c/a)
        bs = BeamSplitter(n, k, ph, rang)
        x!(bs, mat)
    else
        a = mat[n, k]
        rang = atan(abs(c) , abs(a))
        ph = abs2(a) <=0 ? 0 : angle(-c/a)
        bs = BeamSplitter(k, m, ph, rang, true)
        x!(mat, bs)
    end
    return bs
end

"Represents a unitary matrix in Oxford decomposed form."
struct OxfordDecomp
    bsvec::Vector{BeamSplitter}
    diag::Vector{ComplexF64}
end

"""
`OxfordDecomp(mat::Matrix{T<:Number})`

Return an OxfordDecomp object giving the oxford decomposition of
a unitary matrix `mat`.
"""
function OxfordDecomp(mat::Matrix{<:Number})
    resmat = copy(mat)
    bsvec = oxford!(resmat)
    matdiag = diag(resmat)
    OxfordDecomp(bsvec, matdiag)
end

Matrix(od::OxfordDecomp) = compose!(od.bs, Diagonal(od.diag))

#left kills with lines; right with columns
"Decompose the matrix using the Oxford decomposition as in the article
(*cite article here*)
This function transforms a unitary matrix into a unitary diagonal matrix.
It returns an array of the beam spillters used to decompose the matrix."
function oxford!(mat::AbstractMatrix)
    n = size(mat)[1]
    bsplitters = BeamSplitter[]
    left = false
    x = 1
    y = 1
    dir = 1
    while x<n && y <n
        # eliminate entries in the lower triangle in "zigzag" order.
        if left
            bs = elim(mat, n-x+1, y, n-x, left)
        else
            bs = elim(mat, n-x+1, y, y+1, left)
        end
        push!(bsplitters, bs)
        x += dir
        y -= dir
        if x == 0
            x = 1
            dir *= -1
            left = !left
        elseif y == 0
            y = 1
            dir *= -1
            left = !left
        end
    end
    return bsplitters
end

"Compose an interferometer with beam splitters and phase shifters according to
the Oxford prescription. Inverts the Oxford! function if given its output as the
arguments."
function compose!(array::Vector{BeamSplitter}, matrix)
    for j in Iterators.reverse(array)
        j.transposed ? x!(matrix, inverse(j)) : x!(inverse(j), matrix)
    end
    matrix
end

"Compose an array without overwriting it."
function compose(array::Vector{BeamSplitter}, matrix)
    result = copy(matrix)
    compose!(array, result)
end

function compose(array::Vector{BeamSplitter})
    dim = max([max(i.n, i.m) for i in array])
    matrix = Matrix{ComplexF64}(I, dim, dim)
    compose!(array, matrix)
end

"Builds a beam splitter according to the Oxford construction from array of variables.
Intended for use with the optimize function of the Optims.jl package."
function oxfordarray(bsphases, bsmixing, dim)
    nsplits = div(dim*(dim-1), 2)
    vector = Vector{BeamSplitter}(undef, nsplits)
    left = false
    x = 1
    y = 1
    dir = 1
    i = 1
    while x<dim && y <dim
        # eliminate entries in the lower triangle in "zigzag" order.
        if left
            bs = BeamSplitter(dim-x+1, dim-x, bsphases[i], bsmixing[i])
        else
            bs = BeamSplitter(y+1, y,  bsphases[i], bsmixing[i], true)
        end
        vector[i] = bs
        x += dir
        y -= dir
        i += 1
        if x == 0
            x = 1
            dir *= -1
            left = !left
        elseif y == 0
            y = 1
            dir *= -1
            left = !left
        end
    end
    return vector
end

