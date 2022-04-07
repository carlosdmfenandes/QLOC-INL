"Transforms interferometres between the array and
 vector of parameters representations."

include("Reck-Oxford.jl")

"Take a vector of parameter and returns a
dim*dim matrix by oxford composition.
The inverse of argument form."
@views function matrix_form(vec::Vector{<:Real}, dim)
    array_size = (length(vec)-dim)÷2 #determines the number of BeamSplitters
    phvec = vec[1:dim] #Interpret slice as free phase angles
    mixes = vec[dim+1:dim+array_size] #Interpret slice as BS mixing angle
    phases = vec[dim+1+array_size:dim+2*array_size] #Interpret slice as BS phase angle
    bsvec = oxfordarray(phases, mixes, dim)
    dmat= Matrix( Diagonal(cis.(phvec)) )
    compose(bsvec, dmat)
end

"Deduce dimension from the intended matrix form the length of the argument vector."
function matrix_form(vec)
    n = isqrt(length(vec))
    matrix_form(vec, n)
end

function argument_form(matrix)
    "Take a matrix and obtain the vector of parameter by oxford decomposition.
    The inverse of matrix form."
    dim = size(matrix)[1] #declares a vector with n^2 parameter of a unitary matrix
    vec = Vector{Float64}(undef, dim*dim)
    fmat = copy(matrix)
    beams = oxford!(fmat)
    beamn = length(beams)
    fdiag = angle.(diag(fmat))
    bphases = [i.phase for i in beams]
    bmixes = [i.rangle for i in beams]
    vec[1:dim] = fdiag
    vec[dim+1:dim+beamn] = bmixes
    vec[dim+beamn+1:dim+2*beamn] = bphases
    return(vec)
end

#reduced_form(function, n, args...) = x -> function(matrix_form(x,n), args...)
