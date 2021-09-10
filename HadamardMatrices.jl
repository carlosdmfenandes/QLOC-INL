function sylvester(exp)
    base = [1]
    mat = [1 1; 1 -1]
    for i in 1:exp
        base = kron(base,mat)
    end
    return base
end

function fourier_matrix(n)
    matrix = [exp(2*pi*im*i*j/n) for i=0:n-1, j=0:n-1]
    return matrix
end

