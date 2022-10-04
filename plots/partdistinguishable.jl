#=
This file is under the GNU GPL version 2 License.
Copyright 2022 carlosfernandes <carlosfernandes@protonmail.com>
=#
include("../csvread.jl")
include("../partdistinguishable.jl")
using Optim
using DataFrames
using Plots
using Alert #Used to send notifications after terminating calculations.
pgfplotsx()

"Lousy, ad-hoc gram matrix constructor."
function smatt(x, y)
    svec = [1 0 0; x sqrt(1-x^2) 0; y 0 sqrt(1-y^2)]
    svec*adjoint(svec)
end

"Lousy, ad-hoc gram matrix constructor."
function fourmatt(alph, bet, gam, thet)
    svec = [1 0 0;
            cos(alph) sin(alph) 0;
            cos(gam) cis(thet)*sin(beth)*sin(gam) sin(gam)*cos(beth)]
    svec*adjoint(svec)
end

"Lousy, ad-hoc gram matrix constructor."
function matt4(alph, bet, gam, thet)
    svec=[
        1 0 0;
        alph sqrt(1-alph^2) 0;
        gam cis(thet)*sqrt(1-beth^2)*sqrt(1-gam^2) sqrt(1-gam^2)*beth
        ]
    svec*adjoint(svec)
end

function smatt4(x)
    mat = [j==1 ? x : zero(x) for i=1:4, j=1:4] + I*sqrt(1-x^2)
    mat[1,1] = 1
    return mat*adjoint(mat)
end

"""
Define a function that return the no-bunching probability given a 2d vector of
overlaps.
"""
function partopt(umatrix)
    function optf(vec::Vector{Float64})
        partdistinguishable(umatrix, smatt(vec[1],vec[2]))
    end
    return optf
end

function gap(matrix, e=eps())
    upper = [1-e, 1-e]
    lower = [e, e]
    init = [rand(), rand()]
    f = partopt(matrix)
    opt = optimize(f, lower, upper, init, Fminbox(GradientDescent()))
    tup = (f(lower), opt.minimum, f(upper))
    if tup[3] < tup[1]
        delta = max(tup...) - min(tup...)
    else
        delta = min(tup...) - max(tup...)
    end
    (opt.minimizer, delta)
end

function gap4(matrix, tol=eps(), e=eps())
    upper = [1-tol, 1-tol, 1-tol, 2*pi]
    lower = [tol, tol, tol, 0]
    rands = rand(Float64, 4)
    init = lower .+ (upper .- lower).*rands
    f = partopt(matrix)
    opt = optimize(f, lower, upper, init, Fminbox(GradientDescent()))
    tup = (f(lower), opt.minimum, f(upper))
    if tup[3] < tup[1]
        delta = max(tup...) - min(tup...)
    else
        delta = min(tup...) - max(tup...)
    end
    (Tuple(opt.minimizer), delta)
end

function table(mats, tol=eps())
    header = ["xmin", "ymin", "gap"]
    df = DataFrame([Float64[] for i=1:3], header)
    table!(df::DataFrame, mats, tol)
end

function table!(df::DataFrame, mats, tol=eps())
    for m in mats
        g = gap(m, tol)
        row = [g[1][1], g[1][2], g[2]]
        push!(df, row)
    end
    df
end

function table4!(df::DataFrame, mats, tol=eps())
    for m in mats
        g = gap4(m, tol)
        row = [g[1][1], g[1][2], g[1][3], g[1][4], g[2]]
        push!(df, row)
    end
    df
end

function pathgen(dir, range)
    ndir = normpath(dir)
    gen = ((joinpath(ndir, "R$i.csv"), joinpath(ndir, "I$i.csv"))
            for i in range)
end

function checkpath(path)
    verify = ispath(path)
    warnstring = """Cannot construct matrix in $path.\n Necessary files not found."""
    if !verify
        @warn warnstring
    end
    verify
end

matgen(paths) = (matread(i[1], i[2]) for i in paths if all(checkpath, i))

function dirmats(dir, range=0:999)
    matarray = Matrix{ComplexF64}[]
    ndir = normpath(dir)
    for i in range
        rpath = joinpath(ndir, "R$i.csv")
        ipath = joinpath(ndir, "I$i.csv")
        if isfile(rpath) && isfile(ipath)
            mat = matread(rpath, ipath)
            push!(matarray, mat)
        else
            @warn """
                  Cannot construct matrix $i in $dir.
                  Necessary files not found.
                  """
        end
    end
    matarray
end

function datatable(dirs, tol)
    header = ["xmin", "ymin", "gap"]
    df = DataFrame([Float64[] for i=1:3], header)
    for d in dirs
        mats = pathgen(d, 0:999) |> matgen
        table!(df, mats, tol)
    end
    df
end

function datatable4(dirs, tol)
    header = ["amin", "bmin", "gmin", "thmin", "gap"]
    df = DataFrame([Float64[] for i=1:5], header)
    for d in dirs
        mats = pathgen(d, 0:999) |> matgen
        table4!(df, mats, tol)
    end
    df
end
#=
  -0.74411-0.169724im   0.313845+0.354604im   0.437846+0.0392701im
  0.170136-0.424952im  0.0873588-0.625108im   0.587493-0.216632im
 -0.393892-0.22981im    0.284193-0.54461im   -0.551018+0.333243im

  0.467617+0.343499im   -0.483362-0.194312im   0.619199+0.092406im
 -0.614384-0.0567497im  -0.638358+0.39332im   0.0873792+0.222428im
 -0.247409+0.470577im   -0.269987-0.305793im  -0.278235-0.688135im
=#
