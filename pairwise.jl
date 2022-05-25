#=
 * Copyright 2022 carlosfernandes <carlosfernandes@protonmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
=#

"""
Take the next term of the sum according to Glynn's Formula.
Could it be more elegantly implemented with Julia's iteration
interface?
"""
function  iter!(column, icount, matrix)
    length = size(matrix)[1]
    @inbounds if icount != 0
        pos, sign = grayBitToFlip(icount)
        if sign == 1
            @simd for j in 1:length
                column[j] -= matrix[j, pos]
            end
        else
            @simd for j in 1:length
                column[j] += matrix[j, pos]
            end
        end
    end
    prd = one(eltype(column))
    @inbounds @simd for j in 1:length
        prd *= column[j]
    end
    term = partosign(icount)*prd
    icount += 1
    return (term, icount)
end

"""
Do the pairwise summation.
With some refactoring it could be generalized to a lazy pairwise
summation over any iterator.
"""
function pairwise_permantent(matrix)
    length::UInt = size(matrix)[1]
    nterms::UInt = 2^(length-1)
    storage = Vector{eltype(matrix)}(undef, length)
    addr = 0 #points to the largest occupied
    column = sum(matrix, dims=2)
    column ./= 2
    count = 0 #counts the number of iterations passed
    @inbounds for i in 1:((nterms+1)÷2)
        addr += 1
        term, count = iter!(column, count, matrix)
        storage[addr] = term
        if count < nterms
            term, count = iter!(column, count, matrix)
            storage[addr] += term
        end
        for j in 1:trailing_zeros(i)
            addr -= 1
            storage[addr] += storage[addr + 1]
        end
    end
    @inbounds while addr>1
        addr -= 1
        storage[addr] += storage[addr + 1]
    endas
    2*storage[1]
    end
end

function pairwise(terms_iter)
    T::Type=Float64
    len::UInt = length(terms_iter)
    storagesize = 8*sizeof(len) - leading_zeros(len)
    storage = Vector{T}(undef, storagesize)
    addr = 0
#    assign_next = true
    @inbounds for (i,t) in enumerate(terms_iter)
        if isodd(i)
            addr += 1
            storage[addr] = t
        else
            storage[addr] += t
            for j in 1:trailing_zeros(i>>1)
                addr -= 1
                storage[addr] += storage[addr + 1]
            end
        end
    end
    @inbounds while addr>1
        addr -= 1
        storage[addr] += storage[addr + 1]
    end
    storage[1]
end

#=
"""
A incomplete implementation of an iterator

struct Glynniter{Tm<:Number, Ti<:Integer}
    matrix::Matrix{Tm}
    i::Ti
    column::Vector{Tm}
end

function Glynniter(mat::AbstractMatrix, i::Integer)
    column = sum(mat, dims=2)
    Glynniter(mat, i, column)
end

function iterate(it::Glynniter)
    mat = it.matrix
    length = size(mat,1)
    lenfac = 2^(length-1)
    column = sum(mat, dims=2)
    column ./= 2
    total = prod(column)
    pos, sign = grayBitToFlip(i)
    if sign == 1
        @simd for j in 1:length
            column[j] -= mat[j, pos]
        end
    else
        @simd for j in 1:length
            column[j] += mat[j, pos]
        end
    end
    prd = one(eltype(column))
    @simd for j in 1:length
        prd *= column[j]
    end
    total += partosign(i)*prd
end
"""
=#
