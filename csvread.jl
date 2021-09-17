using DelimitedFiles

"Defines utilites to process matrices store in csv format."

"Short convinient function to condense our code."
myread(file) = readdlm(file, ',', Float64)

"Import a matrix from two csv files each containg the real and imaginary
 parts."
function matread(source_r, source_i)
    a = myread(source_r)
    b = myread(source_i)
    omat = a + im .* b
end
