using DelimitedFiles

home = homedir()

source_r = "$home/Documents/Quantum_Information/random_matrices/dim_9/Unitary_29R_example2.csv"
source_i = "$home/Documents/Quantum_Information/random_matrices/dim_9/Unitary_29I_example2.csv"

myread(file) = DelimitedFiles.readdlm(file, ',', Float64)

a = myread(source_r)
b = myread(source_i)
omat = a + im.*b

rot(th) = ComplexF64[cos(th) -sin(th); sin(th) cos(th)]
