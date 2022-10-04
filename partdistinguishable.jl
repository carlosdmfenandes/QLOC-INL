#=
This file is under the GNU GPL version 2 License.
Copyright 2022 carlosfernandes <carlosfernandes@protonmail.com>
=#

#include("QLOC.jl")
include("Permanent.jl")
include("HadamardMatrices.jl")
include("permutes.jl")
#using .QLOC
using .UnsafePermutes

function sumbroadcasted(A)
    res, iter = Iterators.peel(A)
    for i in iter
        res = res .+ i
    end
    res
end

function directsum(args...)
    Types = eltype.(args)
    mT = promote_type(Types...)
    sizes = size.(args)
    msize = sumbroadcasted(sizes)
    res = Matrix{mT}(undef, msize)
    ends = (1, 1)
    for (block, sizetup) in zip(args, sizes)
        starts = ends
        ends = starts .+ sizetup
        foreblock = view(res, starts[1]:(ends[1]-1), 1:(starts[2]-1))
        afterblock = view(res, starts[1]:(ends[1]-1), ends[2]:msize[2])
        fill!(foreblock, 0)
        fill!(afterblock, 0)
        res[starts[1]:(ends[1]-1), starts[2]:(ends[2]-1)] = block
    end
    res
end

#=
dft_mat = directsum(fourier_matrix(5), I(2))
bs1 = BeamSplitter(1, 6, 2*pi*rand(), 2*pi*rand())
bs2 = BeamSplitter(2, 7, 2*pi*rand(), 2*pi*rand())
interf = (dft_mat*bs1)*bs2
interf = directsum(interf, interf)

cols = [(j & 2^(i-1) != 0) ? i+7 : i for i in 1:5, j in 0:(2^5-1)]
argvecs = [[(i==6 || i==14 || i in j) for i in 1:14]
            for j in eachcol(cols)]
=#

"""Calculates"""
function immanent_term(matrix::AbstractMatrix, per)
    matdim = size(matrix)[1]
    total = zero(eltype(matrix))
    for (perm, par) in zip(permutations(1:matdim), params)
        prod = one(eltype(matrix))
        @inbounds @simd for j in 1:matdim
            prod *= matrix[perm[j], j]
        end
        total += prod*par
    end
    return total
end

function partdistinguishable(umatrix, smatrix)
    T = eltype(umatrix)
    matdim = size(umatrix)[1]
    adjmat = adjoint(umatrix)
    total = zero(T)
    @inbounds for p in Permutes(collect(1:matdim))
        sview = view(smatrix, p, :)
        par = permanent(sview .* adjmat)
        prod = one(T)
        @inbounds @simd for j in 1:matdim
            prod *= umatrix[p[j], j]
        end
        total += par*prod
    end
    real(total)
end
