"""A lazy pairwise summation
    could be generalised to any lazy iterator
"""
function pairwise(iter, nterms::Int)
    storagesize = 8*sizeof(nterms) - leading_zeros(nterms)
    storage = Vector{Float64}(undef, storagesize)
    addr = 1 #points to the largest occupied

    function get()
        addr += 1
        storage[addr] = iter()
    end

    function add()
        addr -= 1
        storage[addr] += storage[addr + 1]
    end

    function getadd()
        storage[addr] += iter()
    end

    for i in 1:nterms
        get()
        getadd()
        for j in 1:trailing_zeros(i)
            add()
        end
    while addr>1
       add()
    end
    storage[1]
end

pairwise(iter) = pairwise(iter, length(iter))
