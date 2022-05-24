# Part of QLOC.jl package
using LinearAlgebra
using Printf

"""
A structure representing a beamsplitter with a phase shift acting on two
modes `m` and `n` of a system of unspecified length.
This beam splitter is represented by a unitary matrix of the form:

[s exp(im ph) c ;
 c exp(im ph) s ;]

where `s, c = cis(rangle)` and `ph = phase`. If the `transposed` field
is `true` then the structure reprents the transpostion of the matrix above.

Note that `BeamSplitter <: AbstractArray` is `false`.
This is because `AbstractArray` requires implementation of array size,
which is not possible since the system's length is left unspecified;
c.f. `UniformScaling` in the standard `LinearAlgebra` package.
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

#Define some standard array functions on BeamSplitter

"""
Get the matrix element of a `BeamSplitter` object with `inds`. Since
a beam splitter is represented by a matrix, any indices beyond the first
two are ignored.

Indices are implement according to the following rules:

    `bs[bs.m, bs.m] = cos(θ) * cis(φ)`
    `bs[bs.m, bs.n] = -sin(θ)`
    `bs[bs.n, bs.m] = sin(θ) * cis(φ)`
    `bs[bs.n, bs.n] = cos(θ)`
    `bs[x,y] = δ(x,y)` for x or y out of the beamsplitter.
"""
function getindex(bs::BeamSplitter, inds...)
    vinds = inds[1:2]
    ph = cis(bs.phase)
    s::ComplexF64, c::ComplexF64 = sincos(bs.rangle)
    if vinds[1] == bs.m
        if vinds[2] == bs.m
            res = s * ph
        elseif vinds[2] == bs.n
            res = c
        end
    elseif vinds[1] == bs.n
        if vinds[2] == bs.m
            res = c * ph
        elseif vinds[2] == bs.n
            res = -s
        end
    elseif vinds[1] == vinds[2]
        res = one(ComplexF64)
    else
        res = zero(ComplexF64)
    end
    res
end

"""
Construct the 2*2 unitary Matrix representation of a `BeamSplitter`.
Note however that, in general, `bs*A != Matrix(bs)*A`.
"""
function Matrix(bs::BeamSplitter)
    s, c = sincos(bs.rangle);
    ph = cis(bs.phase);
    mat = ComplexF64[ c*ph -s ;
                      s*ph  c ]
    bs.transposed ? transpose(mat) : mat
end
# Define standard operations on a beam splitter.

transpose(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, bs.phase, bs.rangle, !bs.transposed)
inv(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, -bs.phase, bs.rangle, !bs.transposed)
conj(bs::BeamSplitter) =
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
        LinearAlgebra.axpby!(-s, X, c, Y) #Y = u[m, n]X + u[n, n]Y
        LinearAlgebra.axpby!(s * ph, oldY, c * ph, X) #X = u[m, m]X + u[m, n]Y
    end
    return mat
end

"Multiply a matrix by a beamsplitter to the left in place."
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
*(A::BeamSplitter, B::AbstractArray) = x!(A, copy(B))
*(A::AbstractArray, B::BeamSplitter) = x!(copy(A), B)

"
Convert an array of beam splitters to code for the QTikz LaTeX package
that draws a diagram of the array.
"
function btoQTikz(body::Vector{String})
    preamble = [
        raw"\documentclass[tikz]{standalone}",
        raw"\usepackage{tikz}",
        raw"\usetikzlibrary{quantikz}",
    ]
    beginwrapper = [
        raw"\begin{document}",
        raw"\begin{quantikz}[transparent]",
    ]
    endwrapper = [
        raw"\end{quantikz}",
        raw"\end{document}"
    ]
    join(append!(preamble, beginwrapper, body, endwrapper), '\n')
end

"
Return the body of a qtikz block representing an array `bsvec` of
BeamSplitter objects.
"
function qtikzbody(bsvec, header=[])
    nlines = maximum(bs->max(bs.n, bs.m), bsvec)
    mstr = [[""] for i=1:nlines]
    for (l, h) in zip(mstr, header)
        append!(l,h)
    end
    for bs in bsvec
        for (index, array) in pairs(mstr)
            push!(array, qtikzcell(bs, index))
        end
    end
    [join(push!(strarray, "\\\\"), " & ", raw" & \qw ") for strarray in mstr]
end

"Insert a column representing the action of a single beam splitter."
function qtikzcell(bs::BeamSplitter, row_num::Int)
    low, hi = minmax(bs.n, bs.m)
    anglestring = @sprintf("%4.2f", bs.rangle)
    if row_num == low
        str = "\\gate[$(hi-low+1)]{$(anglestring)}"
    elseif low < row_num < hi
        str = raw"\strikethrough"
    else
        str = raw"\qw"
    end
    str
end

function qtikzcell(bs::BeamSplitter, row_num::Int)
    low, hi = minmax(bs.n, bs.m)
    anglestring = @sprintf("%4.2f", bs.rangle)
    phstring = @sprintf("%4.2f", bs.phase)
    if row_num == low
        if bs.transposed
            str = """
                \\gate[$(hi-low+1)]{$(anglestring)} &&
                 \\gate[$(hi-low+1)]{$(anglestring)}
                """
    elseif low < row_num < hi
        str = raw"\strikethrough && \strikethrough"
    else
        str = raw"\qw && \qw"
    end
    str
end

#=
DeadCode
"""
`diagonality(fnorm::Function, mat::Matrix)`

Return the diagonality of `mat` as given by the norm function `fnorm`.
`fnorm` defaults to `norm` if ommited.
"""
diagonality(fnorm::Function, mat::Matrix) = fnorm(diag(mat))/fnorm(mat)
diagonality(mat::Matrix) = diagonality(norm, mat)

"
Represent a phase shifter acting on a single optical mode `m`,
shifting by a given `phase`.
"
struct PhaseShifter
    phase::Float64
    m::Int
end

"
Represent a simple BeamSplitter mode `m`,
shifting by a given `phase`.
"
struct SimpleBS
    m::Int
    n::Int
    rangle::Float64
end
=#

