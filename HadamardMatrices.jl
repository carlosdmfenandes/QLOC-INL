function sylvester_matrix(exp)
    base = [1]
    mat = [1 1; 1 -1]
    for i in 1:exp
        base = kron(base, mat)
    end
    return base
end

"""
Generates the unnormalized `n*n` DFT matrix F_{jk}=exp(2πi(j-1)*(k-1)/n)
Each root of unity is calculated only once and is copied to the remaining
elements, guaranteeing exact symmetry and equal precision fo all matrix
elements.
"""
function fourier_matrix(n::Integer)
    matrix = Matrix{ComplexF64}(undef, n, n)
    fill!(view(matrix, :, 1), one(ComplexF64))
    if n >= 2
        for i=1:n
            matrix[i, 2] = cispi(2*(i-1)/n)
        end
        for j in 3:n, i in 1:n
            matrix[i,j] = matrix[(i-1)*(j-1)%n+1, 2]
        end
    end
    return matrix
end

function simple_fourier_matrix(n)
    matrix = [cispi(2*i*j/n) for i=0:n-1, j=0:n-1]
    return matrix
end
