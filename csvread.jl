using DelimitedFiles

"Defines utilites to process matrices store in csv format."

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
