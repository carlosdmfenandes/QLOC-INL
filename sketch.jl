function bitoddity(n)
    s = n
    odd = false
    while s != 0
        z = !iszero(s & 1)
        odd = xor(odd, z)
        s = s >>> 1
    end
    return odd
end
randoddity() = bitoddity(rand(Int) & rand(Int))
randoddity(n) = [randoddity() for i in 1:n]
