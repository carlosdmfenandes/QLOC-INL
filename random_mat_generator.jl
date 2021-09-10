using RandomMatrices
using DelimitedFiles

DIR = "$(homedir())/Documents/Quantum_Information/random_matrices/juliarands"
NUM = 100

function randHaarwrite(dim, range=1:NUM, dir=DIR, split=1000)
    writedir = "$dir/dim$dim"
    if !isdir(writedir)
        mkdir(writedir)
    end
    for i in range
        if i%split == 1
	    writedir = "$writedir/$(div(i,split))"
            if !isdir(writedir)
                mkdir(writedir)
            end
        end
        mat = rand(Haar(2), dim)
        rmat = real.(mat)
        imat = imag.(mat)
	print(writedir)
        rfile = "$writedir/$(i)_R.csv"
        ifile = "$writedir/$(i)_I.csv"
        writedlm(rfile, rmat, ',')
        writedlm(ifile, imat, ',')
    end
end
