struct BinaryIndexIterator
    d::Int64#  dimension
    indices::Array{Int64}# indices of 1's
    init_value::Int64# value that the iterator add the subnumbers to it
end

"""
For some a in {0,1}^d, that have 1's in the places "indices", iterate over all substes 
of indices in "indices" and return the number of it's binary representation
"""
function BinaryItr(d::Int64, indices::Array{Int64}, init_value::Int64 = 0)
    @assert 0 < d <= 62

    BinaryIndexIterator(d, indices, init_value)
end

"""
Now the indices of the 1's are the number itself with those indices
"""
function BinaryItr(d::Int64, indices::Int64, init_value::Int64 = 0)
    @assert 0 < d <= 62

    arr_inx = []
    #one_index = 1

    for i=1:d
        if indices == 0
            break
        end
        if indices%2 == 1
            push!(arr_inx, i)
            #one_index += 1
        end
        indices >>= 1
    end

    BinaryIndexIterator(d, Array{Int64}(arr_inx), init_value)
end

const BinaryIndexState = Int64

eltype(it::BinaryIndexIterator) = Int64
length(it::BinaryIndexIterator) = 2^(length(it.indices))

function start(it::BinaryIndexIterator)
   0
end

function iterate(it::BinaryIndexIterator, state::BinaryIndexState = start(it))

    if state == length(it)
        return nothing
    end

    sub_binary_number = it.init_value

    temp_state = state

    for i=1:length(it.indices)
        if temp_state == 0
            break
        end
        if temp_state%2 == 1
            sub_binary_number += 1<<(it.indices[i]-1)
        end
        temp_state >>= 1
    end

    state += 1

    return sub_binary_number, state

end