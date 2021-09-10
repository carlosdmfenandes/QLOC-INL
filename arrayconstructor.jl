function a_stackrape(n, m)
    array = Vector{Int}[]
    if m == 1
        return [[n]]
    end
    for i in 0:n
        midarray = a_stackrape(n-i,m-1)
        for j in midarray
            push!(j, i)
            push!(array, j)
        end
    end
    return array
end

function adder!(array, n, m)
    for i in 0:n
        add = copy(array[m])
        push!(add, i)
        push!(array, i)
    end
end

function a_num(n, m)
    result = Int[]
    interm = [0]
    for i in 1:(n+m-1)
        array = deepcopy(interm)
        for j in array
            a = xor(j, 1 << i-1)
            ons = count_ones(a)
            if ons < n
                push!(interm, a)
            elseif ons == n
                push!(result,a)
            end
        end
    end
    return bitstring.(result)
end

function get_array!(l, array, N, index=1)
    if N == 0
        push!(l, array)
    elseif index == length(array) + 1
        return
    else
        for i in 0:N
            new_array = copy(array)
            new_array[index] += i
            get_array!(l, new_array, N - i, index + 1)
        end
    end
end

function a_matos(n, m)
    res = Vector{Int}[]
    base_array = zeros(Int, m)
    get_array!(res, base_array, n)
    return res
end
