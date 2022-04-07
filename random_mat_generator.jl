import RandomMatrices: Haar
using DelimitedFiles

DIR = "$(homedir())/Documents/Quantum_Information/random_matrices/juliarands"
NUM = 100

function randHaarwrite(dim, range=1:NUM, dir=DIR, dirsize=1000)
    basedir = joinpath("$dir", "dim$dim")
    if !isdir(basedir)
        mkdir(basedir)
    end
    writedir="$basedir"
    for i in range
        if i%dirsize == 1
            writedir = joinpath("$basedir","$(div(i,dirsize)+1)")
            if !isdir(writedir)
                mkdir(writedir)
            end
        end
        mat = rand(Haar(2), dim)
        file = joinpath("$writedir", "$(i).csv")
        writedlm(file, mat, ',')
    end
end
