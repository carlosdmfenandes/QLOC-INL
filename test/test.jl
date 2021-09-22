"""
Benchmark performance of calculationof the permanent.
Demonstrate performance boost for matrices with repeated
lines and columns."""

using BenchmarkTools

include("../Permanent.jl")

test_matrix = rand(ComplexF64, 10, 10)
indices = ones(Int, 20)
indices[11:20] = 1:10
large_matrix = test_matrix[indices, indices]
rinds = ones(Int, 10)
rinds[1] = 11

@btime permanent(large_matrix)
@btime rmultipermanent(test_matrix, Tuple(rinds), Tuple(rinds))
@btime altrmultipermanent(test_matrix, Tuple(rinds), Tuple(rinds))

"""
Benchmark Specs:

HP EliteBook 840 G3
CPU: Intel(R) Core(TM) i7-6500U CPU @ 2.50GHz
RAM: Samsung M471A1G43DB0-CPB DDR4 8GB 2133Mhz 

My results:

First Method    102.870 ms (2 allocations: 432 bytes)
Second Method   885.607 μs (6147 allocations: 1.41 MiB)
Third Method    1.047 ms (4 allocations: 464 bytes)
"""
