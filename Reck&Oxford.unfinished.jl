using LinearAlgebra
using Logging

struct BeamSplitter
    m::Int
    n::Int
    phase::Float64
    rangle::Float64
    transposed::Bool
    function BeamSplitter(m, n, phase, rangle, transposed)
        if m == n
            error("input indices must be different")
        else
            new(m, n, phase, rangle, transposed)
        end
    end
end

BeamSplitter(m, n, phase, rangle) = BeamSplitter(m, n, phase, rangle, false)

transpose(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, bs.phase, bs.rangle, !bs.transposed)
inverse(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, -bs.phase, bs.rangle, !bs.transposed)
conjugate(bs::BeamSplitter) =
    BeamSplitter(bs.m, bs.n, -bs.phase, bs.rangle, bs.transposed)

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

function x!(bs::BeamSplitter, mat::AbstractVecOrMat)
    cb = transpose(bs)
    x!(Transpose(mat), cb)
    return mat
end

⊙(A::AbstractArray, B::AbstractArray) = A * B
⊙(A::BeamSplitter, B::AbstractArray) = x!(A, deepcopy(B))
⊙(A::AbstractArray, B::BeamSplitter) = x!(deepcopy(A), B)

function deone(mat, n=1, m=2, left=true)
    dim = size(mat)[1]
    c = mat[n, m]
    if abs2(c) <= 0.0
        return mat
    elseif left
        a = mat[m, m]
        rang = atan(abs(c) , abs(a))
        ph = abs2(a) <= 0.0 ? 0.0 : angle(-c/a)
        bs = BeamSplitter(m, n, ph, rang)
        display(bs⊙Matrix{ComplexF64}(I, 4, 4))
        x!(bs, mat)
    else
        a = mat[n, n]
        rang = atan(abs(c) , abs(a))
        ph = abs2(a) <=0 ? 0 : -angle(c/a)
        bs = BeamSplitter(m, n, ph, rang, true)
        display(bs⊙Matrix{ComplexF64}(I, 4, 4))
        x!(mat, bs)
    end
end

function elim(mat, n, m, k, left=true)
    dim = size(mat)[1]
    c = mat[n, m]
    if abs2(c) <= 0.0
        return mat
    elseif left
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

function next_fill_curve(sx=1, sy=1, dir=1)
    x = sx + dir
    y = sy - dir
    if x == 0
        x = 1
        dir *= -1
    elseif y == 0
        y = 1
        dir *= -1
    end
    x, y, dir
end

#left kills with lines; right with columns

function oxford(mat::AbstractMatrix)
    n = size(mat)[1]
    bsplitters = BeamSplitter[]
    left = false
    x=  1
    y = 1
    dir = 1
    while x<n && y <n
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

function compose!(array::Vector{BeamSplitter}, matrix)
    for j in Iterators.reverse(array)
        j.transposed ? x!(matrix, inverse(j)) : x!(inverse(j), matrix)
    end
    matrix
end

function compose(array::Vector{BeamSplitter})
    dim = max([max(i.n, i.m) i in array])
    matrix = Matrix{ComplexF64}(I, dim, dim)
    compose!(array::Vector{BeamSplitter}, matrix)
end
