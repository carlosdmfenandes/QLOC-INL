using Statistics
include("Permanent.jl")

"Contains tools for statistical analysis of post-selected nonlinear optical
processes with special attention to merit functions."

"Calculate the success transition amplitude for n input photons
for a given unitary matrix."
function psamplitude(nphotons, matrix)
    dim = size(matrix)[1]
    a = ones(Int, dim)
    a[1] = nphotons
    rmultipermanent(matrix, a, a)/factorial(nphotons)
end

"Generate the matrix relating the polynomial coefficients
to the generated angles."
Lmat(n) = [i^j for i=0:n, j=0:n]

"Calculate the coefficients given a list of angles."
function anglestocoefs(angles)
    n=length(angles)
    inv(Lmat(n))*angles
end

"!A vector of functions? This is strange."
nsamplitudes(n, matrix) = [i -> psamplitude(i, matrix) for i in 0:n]

merit(array, key, args...) = get(MERIT_DICTS, key, nothing)(array, args)

merit_probdiffs(parray) =  sum( abs.( abs.(parray[2:end]) .- abs(parray[1]) ) )#jl arrays start at 1

merit_stddev(parray) = std(abs.(parray), false)

merit_lprob(parray, c=1.) = merit_stddev(parray) + c*(1-mean(parray))^2

"Merit function favouring a small angle"
function merit_centre(parray, c=1., d=1., nangle=2)
    angle = anglestocoefs(parray)[nangle]
    merit_lprob(parray, c) + d*sin(angle)^2
end

"Merit function favouring large angles."
function merit_sides(parray, c=1., d=1., nangle=2)
    angle = anglestocoefs(parray)[nangle]
    merit_lprob(parray, c) + d*cos(angle)^2
end

"Merit function favouring some phiangle."
function merit_angle(parray, c=1., d=1., phiangle,nangle=2)
    angle = anglestocoefs(parray)[nangle]
    merit_lprob(parray, c) + d*sin(angle-phiangle)^2
end

"Merit function that weighs variance exponentially. (Hopefully priotizing it
minimization above other parameters.)"
function merit_exp(parray, a=1., c=1.)
    variance = var(abs.(parray))
    expm1(a*var) + c*(1-mean(parray))^2

"
function merit_value(n, matrix)
    dim = size(matrix)[1]
    p0 = abs(psamplitude(0, matrix))
    p1 = abs(permanent(matrix))
    p2 = abs(psamplitude(2, matrix))
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
"

#reduced_form(functon, n, args..) = x -> function(matrix_form(x,n), args...)

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
