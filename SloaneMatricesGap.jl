"Calculates and prints the permanents of the Matrices in the Sloane's
The reference is http://neilsloane.com/hadamard/
With these permanents the bunching gap is also calculated."

"Calculate 'n!/n^n' which gives the permanent of a classical probability matrix
associated with some interferometer."
classic(n) = prod(i/n for i=1:n-1)

"This Dictionary gives the possible permanents of matrixes of signed ones for
different values of the dimension."
vals = Dict(
    1=>[1],
    2=>[0],
    4=>[8],
    12=>[46080],
    16=>[8028160,50692096,17137664,360448],
    20=>[219414528],
    24=>(2^22*3^2).*[73,567, 585, 969, 1207, 1353, 1609, 2167, 2441,2761, 2889,
                     3145, 3273,3401, 3785, 4041, 4279, 4425, 5449, 5961, 6327,
                     6985, 14007, 18249, 18615, 43191, 116919],
    )

"Transform 'a' into the permanent of a unitary matrix by normalising it by
'sqrt(n^n)'' and then squaring.
This implementation was made in order to
to avoid integer overflow errors when calculating 'n^n'."
bquantum(a,n) = (a/exp(n*log(n)))^2

"Calculate the bunching gap for a Hadamard matrix with permanent 'a'
and dimension 'n'."
gap(a,n) = bquantum(a,n) - classic(n)

"Find the quantum bunching for all possible hadamard matrices"
[map(x -> bquantum(x,i), vals[i]) for i = keys(vals)]

"Find the bunching ratio for the 20D hadamard matrix."
bquantum(219414528,20)/classic(20)

"Find the bunching ratio for all 16D hadamard matrices."
m = [8028160,50692096,17137664,360448]
n = map(x -> bquantum(x,16),m)
n./classic(16)

"Find the bunching ratio for the 12D hadamard matrix."
bquantum(46080,12)/classic(12)

m = (2^22*3^2).*[73,567, 585, 969, 1207, 1353, 1609, 2167, 2441,2761, 2889,3145, 3273,3401, 3785, 4041, 4279, 4425, 5449, 5961, 6327, 6985, 14007, 18249, 18615, 43191, 116919]
n = map(x -> bquantum(x,24),m)
n./classic(24)
