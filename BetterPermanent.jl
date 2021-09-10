using LinearAlgebra

#multPow2(n, m) = n << m
#remPow2(n, m) = n & (multPow2(2, m)-1)
"Implement (-1)^x more efficently."
partosign(x) = 1 - ((x & 1 ) << 1)

"Convert a binary number n to its Gray code representation."
bitToGray(n) = n ⊻ (n>>>1)

"Find the bit flipped in Gray code when going from n-1 to n."
function grayBitToFlip(n::Integer)
	n1 = bitToGray(n)
	n2 = bitToGray(n-1)
	d = n1 ⊻ n2 #find the position of the flip
	sign = 1 - ((iszero(n1 & d)  & 1 ) << 1) #equivalent to (-1)^sign
	j = 0
	while d != 0
		d >>>= 1
		j += 1
	end
	return j, sign
end
"""
Alternative implementation of grayBitToFlip roughly 1.5 times slower

function decompose(number::Integer)
    pos = 0
    x = number
    while true
        x, i = x >>> 1, Bool(x & 1) #equivalent to x,i= divrem(x, 2)
        pos+=1
		if i
			sign = 1 - ((x & 1 ) << 1)
			return pos, sign
		end
	end
end
"""

"Best Optimized implementation of the Permanent"
function permanent(matrix)
	T = eltype(matrix)
    length = size(matrix,1)
    column = zeros(T,length)
    total = zero(T)
    for i in 1:(2^length-1)
		par = partosign(i)
        pos, sign = grayBitToFlip(i)
        LinearAlgebra.axpy!(sign,view(matrix,:,pos), column)
        #column += sign*view(matrix, :, pos)
        total += par*prod(column)
    end
    total = total*partosign(length)
    return total
end

cfac(i, n, sign) = cis(2*pi*im*i/(n+1))*expm1(2*sign*pi*im/(n+1))

function multipermanent(matrix::AbstractMatrix{<:Complex}, fvector, yvector)
	length = size(matrix, 1)
    column = sum(matrix, dims=2) #optimize to major order
    total = prod(column)
    for i in 1:(2^length-1)
        pos, sign = multiGrayBitToFlip(i, fvector)
		fac = cfac(pos, fvector[pos], sign)
        LinearAlgebra.axpy!(fac, view(matrix,:,pos), column)
        #column += sign*view(matrix, :, pos)
        total += par*prod(column .^ yvector)
    end
    total = total/prod(fvector .+ 1)
    return total
end

function multipermanent(matrix::AbstractMatrix{<:Complex}, fvector)
	length = size(matrix,1)
    column = sum(matrix, dims=2) #optimize to major order
    total = prod(column)
	prodfac = prod(fvector .+ 1)
	rootfac = 1
    for i in 1:(prodfac-1)
        pos, sign = multiGrayBitToFlip(i, fvector)
		rootmult = fvector[pos]
		rootfac *= cis(2*pi*im*sign/(rootmult+1))
		fac = cfac(pos, rootmult, sign)
        LinearAlgebra.axpy!(fac,view(matrix,:,pos), column)
        #column += sign*view(matrix, :, pos)
        total += rootfac*prod(column)
    end
    total = total/prodfac
    return total
end

multipermanent(matrix, args...) = multipermanent(ComplexF64.(matrix), args...)

function multiGrayBitToFlip(x, vector=ones(Int, 64))
    for (pos, n) in enumerate(vector)
        if x%(n+1) != 0
            sign = 1 - 2(div(x,n+1)%2)
            return pos, sign
        end
        x = div(x,n+1)
    end
end

mydiv(a::Int, b::Int) = div(a,b)
mydiv(a, b) = a/b

"Permanent by Glynn's formula"
function gpermanent(matrix)
    lenfac = 2^(size(matrix,1)-1)
	column = sum(matrix, dims=2) #optimize to major order
    total = prod(column)
	par = 1
    for i in 1:(lenfac-1)
		par *= -1
        pos, sign = grayBitToFlip(i)
		fac = -2*sign
        LinearAlgebra.axpy!(fac, view(matrix,:,pos), column)
        #column += sign*view(matrix, :, pos)
        total += par*prod(column)
    end
    total = mydiv(total,lenfac)
    return total
end
