
using DelimitedFiles
using Optim

include("matrix_composer.jl")

const HOME = homedir()
const PROJDIR = "$HOME/Documents/Quantum_Information"
const DIM = 4
const FILESDIR = "$PROJDIR/random_matrices/dim_$DIM"

"Generate the matrix relating the polynomial coefficients to the generated angles."
Lmat(n) = [i^j for i=0:n, j=0:n]

"Calculate the coefficients given a list of an"
function anglestocoefs(angles)
    n=length(angles)

"Contains a list of functions and what to name the files
this algorithm generates."
num = "2"
ang = pi / 4

vals_dict = Dict(
    0 => (
        "$PROJDIR/random_matrices/PS$(DIM)stats$num.csv",
        "$PROJDIR/random_matrices/pp$DIM",
        merit_func(DIM),
    ),
    1 => (
        "$PROJDIR/random_matrices/PS$(DIM)stats_dev$num.csv",
        "$PROJDIR/random_matrices/ps$(DIM)_dev",
        merit_func_dev(DIM),
    ),
    2 => (
        "$PROJDIR/random_matrices/PS$(DIM)stats_lprob$num.csv",
        "$PROJDIR/random_matrices/ps$(DIM)_lprob",
        merit_func_lprob(DIM, 0.1),
    ),
    3 => (
        "$PROJDIR/random_matrices/PS$(DIM)stats_centre$num.csv",
        "$PROJDIR/random_matrices/ps$(DIM)_centre",
        merit_func_centre(DIM),
    ),
    4 => (
        "$PROJDIR/random_matrices/PS$(DIM)stats_sides$num.csv",
        "$PROJDIR/random_matrices/ps$(DIM)_sides",
        merit_func_sides(DIM),
    ),
    5 => (
        "$PROJDIR/random_matrices/PS$(DIM)stats_ang$num.csv",
        "$PROJDIR/random_matrices/ps$(DIM)_ang",
        merit_func_ang(DIM, ang, 0.01, 0.1),
    ),
    6 => (
        "$PROJDIR/random_matrices/PS$(DIM)stats_exp$num.csv",
        "$PROJDIR/random_matrices/ps$(DIM)_exp$num",
        merit_func_exp(DIM, 100, 1),
    ),
    7 => (
        "$PROJDIR/random_matrices/PS$(DIM)stats_dummy$num.csv",
        "$PROJDIR/random_matrices/ps$(DIM)_dummy$num",
        x -> abs2(tr(matrix_form(x, DIM))),
    ),
)

FILE, MATDIR, to_optimise = vals_dict[6]

"Short convinient function to condense our code."
myread(file) = readdlm(file, ',', Float64)

"Determine the name of the file to export."
csvname(dim, complex, index) = "Unitary_2$(dim)$(complex)_example$index.csv"

"Import a matrix from two csv files each containg the real and imaginary
 parts."
function matread(rarg, iarg)
    source_r = "$FILESDIR/$rarg"
    source_i = "$FILESDIR/$iarg"
    a = myread(source_r)
    b = myread(source_i)
    omat = a + im .* b
end

function writetofile(stream, index)
    rpath = csvname(DIM, 'R', index)
    ipath = csvname(DIM, 'I', index)
    mat = matread(rpath, ipath)
    res = optimize(to_optimise, argument_form(mat))
    min = Optim.minimizer(res)
    minmat = matrix_form(min, DIM)
    p0m = p0(minmat)
    p1m = p1(minmat)
    p2m = p2(minmat)
    expc = (abs(p0m) + abs(p1m) + abs(p2m)) / 3
    var = ((abs(p0m) - expc)^2 + (abs(p1m) - expc)^2 + (abs(p2m) - expc)^2) / 3
    relative_deviation = sqrt(var) / expc
    nlarg = angle(p2m * p0m / p1m^2) / 2
    prob = abs2(expc)
    writedlm(stream, [index relative_deviation nlarg prob], ',')
end

function writetofile(stream1, stream2, index)
    rpath = csvname(DIM, 'R', index)
    ipath = csvname(DIM, 'I', index)
    mat = matread(rpath, ipath)
    res = optimize(merit_func(DIM), argument_form(mat))
    min = Optim.minimizer(res)
    minmat = matrix_form(min, DIM)
    p0m = p0(minmat)
    p1m = p1(minmat)
    p2m = p2(minmat)
    expc = (abs(p0m) + abs(p1m) + abs(p2m)) / 3
    var = ((abs(p0m) - expc)^2 + (abs(p1m) - expc)^2 + (abs(p2m) - expc)^2) / 3
    relative_deviation = sqrt(var) / expc
    nlarg = angle(p2m * p0m / p1m^2) / 2
    prob = abs2(expc)
    writedlm(stream1, [index relative_deviation nlarg prob], ',')
    writedlm(stream2, minmat, ',')
end

filestream = open(FILE, "w")
try
    for x = 1001:100000
        #       matrix_filepath = "$MATDIR/$x.csv"
        #       matrix_filestream = open(matrix_filepath, "w+")
        #        try
        #            writetofile(filestream, matrix_filestream, x)
        writetofile(filestream, x)
        #       finally
        #            close(matrix_filestream)
        #        end
    end
finally
    close(filestream)
end
print("done.")
