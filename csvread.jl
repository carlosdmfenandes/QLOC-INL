using DelimitedFiles

"Defines utilites to process matrices store in csv format."

"Short convinient function to condense our code."
myread(file) = readdlm(file, ',', Float64)

"Import a matrix from two csv files each containg the real and imaginary
 parts."
function matread(source_r, source_i)
    a::Matrix{Float64}  = myread(source_r)
    b::Matrix{Float64}  = myread(source_i)
    omat = a + im * b
end

function csvtoh5(dir)
    for i in 0:999
        rmat = myread("$dir/R$i.csv")
        imat = myread("$dir/I$i.csv")
        cmat = rmat + im * imat
        h5write("$dir/$i.h5","randHaar",cmat)
    end
end
