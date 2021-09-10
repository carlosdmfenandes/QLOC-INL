condrev(x, cond)= cond ? reverse(x) : x

function lorder(vector)
    #create an empty vector with the same type as the input.
    res = Vector{eltype(vector)}[]
    isempty(vector) && return Vector{eltype(vector)}[[]]
    part = lorder(vector[1:end-1])
    for i in 0:vector[end]
        for j in condrev(part, isodd(i))
            add = copy(j)
            push!(add, i)
            push!(res, add)
        end
    end
    return res
end

function multigray(x, vector=ones(Int, 64))
    for (pos, n) in enumerate(vector)
        if x%(n+1) != 0
            sign = 1 - 2(div(x,n+1)%2)
            return pos, sign
        end
        x = div(x,n+1)
    end
end

function altgen(vector=ones(Int, 64))
    state = zeros(eltype(vector), length(vector))
    res = [copy(state)]
    len = prod(x->x+1, vector)-1
    for i in 1:len
        pos, sign = multigray(i, vector)
        state[pos] += sign
        push!(res, copy(state))
    end
    return res
end
