using DataFrames

"""
get the samples (in dim d) and the domain sizes of
the samles (the i'th  axis of sample has domain size of domain_dim[i])
function histogram_from_sample(samples::AbstractArray, domain_dim::Tuple)
    @assert 0 <= minimum(samples)
    @assert maximum(samples) <= 1
    ones_array = ones(Int64, domain_dim)
    histogram_size = prod(domain_dim)
    h = zeros(Int64, domain_dim)

    function f(row)
        [row[i]*domain_dim[i] for i=1:length(row)]
    end

    new_samples = [f(s)[1] for s in samples]
    for x in round.(Int64, new_samples)

        @inbounds h[x+ones_array] += 1
    end
    collect(1:histogram_size)/histogram_size, h[:]
end
"""

function get_flatten_index(index, dim)
    counter = 0
    pow = 1
    for (i,s) in zip(index,dim)
        counter += i*pow
        pow *= s
    end
    counter
end

function reverse_flatten_index(index, dim)

    new_index = zeros(Int64, length(dim))
    reverse_count = length(dim)

    for (i,s) in enumerate(dim)
        new_index[i] += index % s
        index -= index % s
        index /= s
    end

    new_index

end

"""
get the samples (in dim d) and the domain sizes of
the samles (the i'th  axis of sample has domain size of domain_dim[i])
"""
function histogram_from_sample(samples::AbstractArray, domain_dim::Tuple)
    @assert 0 <= minimum(samples)
    @assert maximum(samples) <= 1
    ones_array = ones(Int64, domain_dim)
    histogram_size = prod(domain_dim)
    h = zeros(Int64, histogram_size)
    rows, cols = size(samples)

    new_samples = similar(samples)

    for i=1:rows
        for j=1:cols
            new_samples[i,j] = samples[i,j]*domain_dim[j]
        end
    end

    for i =1:rows
        x = round.(Int64, new_samples[i, :])
        index = get_flatten_index(x,domain_dim)
        if index == 0
            index = 1
        end
        @inbounds h[index] += 1
    end
    collect(1:histogram_size)/histogram_size, h[:]
end

"""
get the samples (in dim d) and the domain sizes of
the samles (the i'th  axis of sample has domain size of domain_dim[i])
"""
function histogram_from_sample2(samples, num_bins)
    @assert 0 <= minimum(samples)
    @assert maximum(samples) <= 1
    h = zeros(num_bins)
    new_samples = samples*num_bins
    for x in round.(Int64, new_samples)
        @inbounds h[x] += 1
    end
    collect(1:num_bins)/num_bins, h
end

"""
get dataset as a DataFrame and convert it to matrix.
"""
function covertData(data::DataFrames.DataFrame)

    @assert length(size(data)) <= 2
    @assert length(size(data)) >= 1

    if length(size(data)) == 2
        rows, cols = size(data)
    else
        rows = size(data)[1]
        cols = 1
    end

    data_matrix = zeros(rows, cols)

    if length(size(data)) == 2
        for i=1:cols
            a = Array(collect(data[i]))
            data_matrix[1:length(a),i] = a/Float64(maximum(a))
        end
    else
        a = Array(collect(data))
        data_matrix[1:length(a),1] = (a/Float64(maximum(a)))
    end

    data_matrix
end


"""
get a list of values, and return how many times the each value is found.
return: dictionary{value : counter}
"""
function countmemb(itr)
    d = Dict{eltype(itr), Int64}()
    for val in itr
        if isa(val, Number) && isnan(val)
            continue
        end
        d[val] = get!(d, val, 0) + 1
    end
    return d
end

function getDomainDim(data::AbstractArray)

    @assert length(size(data)) <= 2
    @assert length(size(data)) >= 1

    if length(size(data)) == 2
        rows, cols = size(data)
    else
        rows = size(data)[1]
        cols = 1
    end

    domain_dim = zeros(Int64, cols)

    if cols == 1
        domain_dim[1] = length(countmemb(data[:]))
    else
        for i=1:cols
            domain_dim[i] = length(countmemb(data[:,i]))
        end
    end

    Tuple(domain_dim)

end

function getDomainDim(data::DataFrames.DataFrame)

    @assert length(size(data)) <= 2
    @assert length(size(data)) >= 1

    if length(size(data)) == 2
        rows, cols = size(data)
    else
        rows = size(data)[1]
        cols = 1
    end

    domain_dim = zeros(Int64, cols)

    for i=1:cols
        domain_dim[i] = length(countmemb(collect(data[i])))
    end

    Tuple(domain_dim)

end

function getDomainDim2(data::DataFrames.DataFrame)

    @assert length(size(data)) <= 2
    @assert length(size(data)) >= 1

    if length(size(data)) == 2
        rows, cols = size(data)
    else
        rows = size(data)[1]
        cols = 1
    end

    domain_dim = zeros(Int64, cols)

    for i=1:cols
        domain_dim[i] = maximum(Array(collect(data[i])))+1
    end

    Tuple(domain_dim)

end

function readContingencyTable(data::DataFrames.DataFrame, dimension::Int64)

    @assert length(size(data)) == 2

    n,temp = size(data)

    @assert temp == 2
    @assert n > 0

    hist = zeros(Float64, 2^dimension)

    for i in 1:n
        dig = digits(data[[:1]][1][i], base = 10)
        index = sum([dig[k]*2^(k-1) for k=1:length(dig)])
        hist[index] += data[[:2]][1][i]
    end

    hist
end

function getMapping(nomain_size::Int64)

    

end