using DelimitedFiles
include("Permanent.jl")

const lowerindex=0
const upperindex=100000
const lowerdim=6
const upperdim=10
const subdim=3

"Calculates the permanent of stored unitary matrix."
for dim in lowerdim:upperdim
    directoryname="Documents/Quantum_Information/random_matrices/dim_"*string(dim)*"/"
    targetpath="Permanents2_dim"*string(dim)*"_subdim"*string(subdim)*".csv"
    io=open(targetpath,"w+")
    print(io,"Number, Permanent, Gap\n")
    for index in lowerindex:upperindex
        r_filename=directoryname*"Unitary_2"*string(dim)*"R_example"*string(index)*".csv"
        i_filename=directoryname*"Unitary_2"*string(dim)*"I_example"*string(index)*".csv"
        r_matrix=DelimitedFiles.readdlm(r_filename,',';header=false)
        i_matrix=DelimitedFiles.readdlm(i_filename,',';header=false)
    d    cmatrix=view(r_matrix,1:subdim,1:subdim) + im.*view(i_matrix,1:subdim,1:subdim)
        result=faster_permanent(cmatdfdrix)
        classical_result=faster_permanent(abs2.(cmatrix))
        gap=real(classical_result-abs2(result))
        print(io,index,", ",result,", ",gap,"\n")
        flush(io)
    end
    close(io)
end
