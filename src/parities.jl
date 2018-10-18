"""
    Parities

Implementation of parity queries using Fast Hadamard Walsh Transform and Gosper
iteration.
"""
struct Parities <: Queries
    dimension::Int
    order::Int
    idx::Array{Int, 1}
end

"""
    hadamard_basis_vector(index, dimension)

Returns hadamard basis vector given by `index`.
"""
function hadamard_basis_vector(index::Int, dimension::Int)
    hadamard = zeros(Float64, 1 << dimension)
    hadamard[1] = 1.0
    for i = 0:dimension-1
        sign = (index & (1 << i)) > 0 ? -1.0 : 1.0
        @simd for j = 1:(1 << i)
            @inbounds hadamard[j + (1 << i)] = sign * hadamard[j]
        end
    end
    hadamard
end

function Parities(dimension, order)
    idx = [1]
    for r = 1:order
        for s = gosper(dimension, r)
            push!(idx, s+1)
        end
    end
    Parities(dimension, order, idx)
end

function get(queries::Parities, i::Int)
    HistogramQuery(hadamard_basis_vector(queries.idx[i]-1, queries.dimension))
end

function get_parities(queries::Parities, i::Int)
    HistogramQuery(hadamard_basis_vector(queries.idx[i]-1, queries.dimension))
end

function evaluate(queries::Parities, h::Histogram)
    2^queries.dimension * fwht_natural(h.weights)[queries.idx]
end

function initialize(queries::Parities, data::Tabular, ps::MWParameters)
    initialize(queries, Histogram(data), ps)
end


struct FactorParity <: FactorHistogramQuery
    attributes::Array{Int, 1}
end

"""
    FactorParities

Implementation of parity queries over factored histograms.
"""
struct FactorParities <: FactorHistogramQueries
    queries::Array{FactorParity, 1}
end

attributes(q::FactorParity) = q.attributes
default_value(q::FactorParity) = 0.0
get(qs::FactorParities, i::QueryIndex) = qs.queries[i]

function restrict(q::FactorParity, attributes::Array{Int, 1})
    idx = 0
    d = length(attributes)
    for a in q.attributes
        i = findfirst(x -> x==a, attributes)#findfirst(attributes, a)
        idx += 2^(d-i)
    end
    HistogramQuery(hadamard_basis_vector(idx, d))
end

function FactorParities(dimension, order)
    qs = FactorParity[]
    for k = 1:order
        append!(qs, [FactorParity(s) for s = subsets(collect(1:dimension), k)])
    end
    FactorParities(qs)
end

function parity(x::Array{Float64, 1}, attributes::Array{Int64, 1})
    (-1.0)^sum(x[attributes])
end

function evaluate(q::FactorParity, table::Tabular)
    s = 0.0
    n = size(table.data)[2]
    a = attributes(q)
    @simd for i = 1:n
        @inbounds s += parity(table.data[:, i], a)
    end
    s/n
end

function evaluate(qs::FactorParities, table::Tabular)
    evals = zeros(Float64, length(qs.queries))
    @simd for i in 1:length(qs.queries)
        @inbounds evals[i] = evaluate(qs.queries[i], table)
    end
    evals
end

function evaluate(qs::FactorParities, d::FactorHistogram)
    [evaluate(q, d) for q in qs.queries]
end
