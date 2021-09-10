"Transforms interferometres between the array and
 vector of parametert representations."

include("Reck&Oxford.jl")
include("Permanent.jl")

function matrix_form(vec::Vector{Float64}, dim)
    "Take a vector of parameter and returns a
    dim*dim matrix by oxford composition.
    The inverse of argument form."
    array_size = (length(vec)-dim)÷2 #determines the number of BeamSplitters
    phvec = vec[1:dim] #Interpret slice as free phase angles
    mixes = vec[dim+1:dim+array_size] #Interpret slice as BS mixing angle
    phases = vec[dim+1+array_size:dim+2*array_size] #Interpret slice as BS phase angle
    bsvec = oxfordarray(phases, mixes, dim)
    dmat= Matrix( Diagonal(cis.(phvec)) )
    compose(bsvec, dmat)
end

function matrix_form(vec::Vector{Float64})
    "Deduce dimension fo the inteded matrix form the length of the argument
    vector."
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

function p0(matrix)
    dim = size(matrix)[1]
    a = ones(Int, dim)
    a[1] = 0
    rmultipermanent(matrix, a, a)
end

function p2(matrix)
    dim = size(matrix)[1]
    a = ones(Int, dim)
    a[1] = 2
    rmultipermanent(matrix, a, a)/2
end

p1(matrix) = permanent(matrix)

function merit_value(matrix)
    dim = size(matrix)[1]
    a1 = ones(Int,dim)
    a0 = copy(a1)
    a0[1] = 0
    a2 = copy(a1)
    a2[1] = 2
    p0 = abs(rmultipermanent(matrix, a0, a0))
    p1 = abs(permanent(matrix))
    p2 = abs(rmultipermanent(matrix, a2 ,a2))/2
    abs(p0-p1) + abs(p0-p2)
end

function merit_value_dev(matrix)
    dim = size(matrix)[1]
    a1 = ones(Int,dim)
    a0 = copy(a1)
    a0[1] = 0
    a2 = copy(a1)
    a2[1] = 2
    p0 = abs(rmultipermanent(matrix, a0, a0))
    p1 = abs(permanent(matrix))
    p2 = abs(rmultipermanent(matrix, a2 ,a2))/2
    mean = (p0 + p1 + p2)/3
    abs(p0-mean) + abs(p1-mean) + abs(p2 - mean)
end

function merit_value_lprob(matrix, c=1.)
    pr0 = abs(p0(matrix))
    pr1 = abs(p1(matrix))
    pr2 = abs(p2(matrix))
    mean = (pr0 + pr1 + pr2)/3
    variance = ((pr0-mean)^2+(pr1-mean)^2+(pr2-mean)^2)/3
    variance + c*(1-mean)
end

function merit_value_centre(matrix, c=1., d=1.)
    p0m = p0(matrix)
    p1m = p1(matrix)
    p2m = p2(matrix)
    pr0 = abs(p0m)
    pr1 = abs(p1m)
    pr2 = abs(p2m)
    mean = (pr0 + pr1 + pr2)/3
    variance =  ((pr0-mean)^2+(pr1-mean)^2+(pr2-mean)^2)/3
    nlarg = angle(p2m*p0m/p1m^2)/2
    variance + c*(1-mean)^2 + d*sin(nlarg)^2
end

function merit_value_sides(matrix, c=1., d=1.)
    p0m = p0(matrix)
    p1m = p1(matrix)
    p2m = p2(matrix)
    pr0 = abs(p0m)
    pr1 = abs(p1m)
    pr2 = abs(p2m)
    mean = (pr0 + pr1 + pr2)/3
    variance =  ((pr0-mean)^2+(pr1-mean)^2+(pr2-mean)^2)/3
    nlarg = angle(p2m*p0m/p1m^2)/2
    variance + c*(1-mean)^2 + d*cos(nlarg)^2
end

function merit_value_ang(matrix, ang, c=1, d=1)
    p0m = p0(matrix)
    p1m = p1(matrix)
    p2m = p2(matrix)
    pr0 = abs(p0m)
    pr1 = abs(p1m)
    pr2 = abs(p2m)
    mean = (pr0 + pr1 + pr2)/3
    variance =  ((pr0-mean)^2+(pr1-mean)^2+(pr2-mean)^2)/3
    nlarg = angle(p2m*p0m/p1m^2)/2
    variance + c*(1-mean)^2 +d*sin(nlarg-ang)^2
end

function merit_value_exp(matrix, a=1, c=0)
    p0m = p0(matrix)
    p1m = p1(matrix)
    p2m = p2(matrix)
    pr0 = abs(p0m)
    pr1 = abs(p1m)
    pr2 = abs(p2m)
    mean = (pr0 + pr1 + pr2)/3
    variance =  ((pr0-mean)^2+(pr1-mean)^2+(pr2-mean)^2)/3
    expm1(a*variance) + c*(1-mean)
end

reduced_form(functon, n, args..) = x -> function(matrix_form(x,n), args...)

#merit_func(n) = x -> merit_value(matrix_form(x,n))
#merit_func_dev(n) = x -> merit_value_dev(matrix_form(x,n))
#merit_func_lprob(n, c=1) = x -> merit_value_lprob(matrix_form(x,n), c)
#merit_func_centre(n, c=0, d=1) = x -> merit_value_centre(matrix_form(x,n), c, d)
#merit_func_sides(n, c=0, d=1) = x -> merit_value_sides(matrix_form(x,n), c, d)
#merit_func_ang(n, ang, c=0, d=1) = x -> merit_value_ang(matrix_form(x,n), ang, c, d)
#merit_func_exp(n, a=1, c=0) = x -> merit_value_exp(matrix_form(x,n), a, c)

#DIM = 4
#INIT_VAL = ComplexF64[1 1 1 1; 1 -1 1 -1; 1 1 -1 -1; 1 -1 -1 1]/2

#res = optimize(merit_func(DIM), argument_form(INIT_VAL))
#show(res)
#min = Optim.minimizer(res)
#minmat = matrix_form(min, DIM)
#p0m = p0(minmat)
#p1m = p1(minmat)
#p2m = p2(minmat)
#print("p0:$p0m \n p1:$p1m \n p2:$p2m \n")
#display(minmat)
