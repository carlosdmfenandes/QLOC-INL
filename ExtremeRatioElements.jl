using DelimitedFiles
using LinearAlgebra
using StatsPlots

q3= findmin(
        readdlm("/home/carlosfernandes/Documents/Quantum_Information/random_matrices/expectedN_3.csv",',')[:,4])
q4= findmin(
        readdlm("/home/carlosfernandes/Documents/Quantum_Information/random_matrices/expectedN_4.csv",',')[:,4])
q5= findmin(
        readdlm("/home/carlosfernandes/Documents/Quantum_Information/random_matrices/expectedN_5.csv",',')[:,4])

println(q3,q4,q5)

file3r="/home/carlosfernandes/Documents/Quantum_Information/random_matrices/dim_3/Unitary_23R_example65118.csv"
file3i="/home/carlosfernandes/Documents/Quantum_Information/random_matrices/dim_3/Unitary_23I_example65118.csv"
file4r="/home/carlosfernandes/Documents/Quantum_Information/random_matrices/dim_4/Unitary_24R_example43196.csv"
file4i="/home/carlosfernandes/Documents/Quantum_Information/random_matrices/dim_4/Unitary_24I_example43196.csv"
file5r="/home/carlosfernandes/Documents/Quantum_Information/random_matrices/dim_5/Unitary_25R_example99788.csv"
file5i="/home/carlosfernandes/Documents/Quantum_Information/random_matrices/dim_5/Unitary_25I_example99788.csv"

m3= readdlm(file3r,',') + im*readdlm(file3i,',')    
m4= readdlm(file4r,',') + im*readdlm(file4i,',') 
m5= readdlm(file5r,',') + im*readdlm(file5i,',')

vals3 = readdlm("/home/carlosfernandes/Documents/Quantum_Information/random_matrices/expectedN_3.csv",',')[:,4]
vals4 = readdlm("/home/carlosfernandes/Documents/Quantum_Information/random_matrices/expectedN_4.csv",',')[:,4]
vals5 = readdlm("/home/carlosfernandes/Documents/Quantum_Information/random_matrices/expectedN_5.csv",',')[:,4]

h3 = histogram(vals3)
h4 = histogram(vals4)
h5 = histogram(vals5)
savefig(h3,"~/Desktop/hist3.png")
savefig(h4,"~/Desktop/hist4.png")
savefig(h5,"~/Desktop/hist5.png")

savefig(,"~/Desktop/hist.png")
