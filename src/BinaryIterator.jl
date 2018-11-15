struct BinaryIndexIterator
    d::Int64#  dimension
    indices::Array{Int64}# indices of 1's
end

const BinaryIndexState = Int64

eltype(it::BinaryIndexIterator) = Int64
length(it::BinaryIndexIterator) = 2^(length(it.indices))

function start(it::BinaryIndexIterator)
   0
end

function iterate(it::BinaryIndexIterator, state::BinaryIndexState = start(it))

    if state == length(it.indices)
        return nothing
    end

    sub_binary_number = 0

    for i=1:length(it.indices)
        if state == 0
            break
        end
        if state%2 == 1
            sub_binary_number += (1<<it.indices[i])
        end
    end

    state += 1

    return sub_binary_number

end