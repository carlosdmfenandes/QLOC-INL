module UnsafePermutes

import Base: iterate, eltype, length

export Permutes, Permutation, PairTransposition, PermTuple
export NearTransposition

"Abstract supertype of all types representing permutations."
abstract type Permutation end

"""
Represents a transposition of ith element of a collection with the
(i+1)th.
"""
struct NearTransposition <: Permutation
    i::Integer
end

"""
Represents the transposition of ith element of a collection with the
jth element.
"""
struct PairTransposition <: Permutation
    i::Integer
    j::Integer
end

"""
Represents the transposition of an array with N elements.
"""
struct PermTuple{N, Ti<:Integer} <: Permutation
    perm::NTuple{N, Ti}
end

"""
Find the position of the transposition that generates the (i+1)th
element of the set of permutations of N elements.
The transposition exchanges the i and i+1 elements.
"""
function transposition(i, N)
    ret = zero(N)
    d=N
    q, r = divrem(i, N)
    p = isodd(q)
    while iszero(r)
        d -= 1
        p && (ret += 1)
        q, r = divrem(q, d)
        p = isodd(q)
    end
    ret += p ? r : (d - r)
end

"""
An unsafe Iterator wrapper that wraps around an Array and gives
permutations of it when iterated over.
For the sake of performance the iterator stores a single array in
memory and mutates over each iterattion, returning a reference to the
stored array. This means that if that doing any
mutating operation on the return value also mutates the stored array.
Such a mutation will also necessarily affect the contents of the array
in all following iterations, "contaminating" the result over all the
loop.
If the size is shorter than the array length, the iteration is over
the first size elements of the array. If size is longer than the array
length, the iteration will try to access an out-of-bounds element of
the array, leading to a BoundsError if bound checks are enabled or
undefined behaviour otherwise.
"""
struct Permutes{A<:AbstractArray}
    array::A
    size::Int
end

Permutes(a::AbstractArray) = Permutes(a, convert(Int, length(a)) )
eltype(p::Permutes{A}) where {A} = A
length(p::Permutes) = factorial(p.size)

@inline iterate(p::Permutes) = (p.array, 0)
@inline function iterate(p::Permutes, state)
    newstate=state+1
    newstate>=factorial(p.size) && (return nothing)
    j = transposition(newstate, p.size)
    @inbounds begin
        a = p.array[j]
        p.array[j] = p.array[j+1]
        p.array[j+1] = a
    end
    return p.array, newstate
end

end

module SafePermutes

import Base: iterate, eltype, length

export Permutes, Permutation, PairTransposition, PermTuple
export NearTransposition

"Abstract supertype of all types representing permutations."
abstract type Permutation end

"""
Represents a transposition of ith element of a collection with the
(i+1)th.
"""
struct NearTransposition <: Permutation
    i::Integer
end

"""
Represents the transposition of ith element of a collection with the
jth element.
"""
struct PairTransposition <: Permutation
    i::Integer
    j::Integer
end

"""
Represents the transposition of an array with N elements.
"""
struct PermTuple{N, Ti<:Integer} <: Permutation
    perm::NTuple{N, Ti}
end

"""
Find the position of the transposition that generates the (i+1)th
element of the set of permutations of N elements.
The transposition exchanges the i and i+1 elements.
"""
function transposition(i, N)
    ret = zero(N)
    d = N
    q, r = divrem(i, N)
    p = isodd(q)
    while iszero(r)
        d -= 1
        p && (ret += 1)
        q, r = divrem(q, d)
        p = isodd(q)
    end
    ret += p ? r : (d - r)
end

struct Permutes{A<:AbstractArray, N<:Integer}
    array::A
end

Permutes(a::AbstractArray) = Permutes(a, convert(Int, length(a)) )
eltype(p::Permutes{A, N}) where {A, N} = A
length(p::Permutes{A, N}) where {A, N} = factorial(N)

@inline function iterate(p::Permutes{A, N}) where {A, N}
    E = eltype(p.array)
    ntup = NTuple{N, E}(p.array)
    (ntup, 0)
end
@inline function iterate(p::Permutes{A, N}, state) where {A,N}
    newstate=state+1
    newstate>=length(p) && (return nothing)
    j = transposition(newstate, N)
    @inbounds begin
        a = p.array[j]
        p.array[j] = p.array[j+1]
        p.array[j+1] = a
    end
    E = eltype(p.array)
    ntup = NTuple{N, E}(p.array)
    return ntup, newstate
end

end
