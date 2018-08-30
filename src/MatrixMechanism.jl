#include("PrivateMultiplicativeWeights.jl")

#using PrivateMultiplicativeWeights

"""
return the lower bound for the "Query andWorkload Error"
    according to the article:
    "Ecient Batch Query Answering Under Dierential Privacy"
"""
function getLowerBound(queries::HistogramQueries, eps, delta)
    m,n = size(queries.queries)
    P_eps_delta = (2*log2(2/delta))/(eps^2)
    singular_values = svdfact(queries.queries)[:S]
    ss_svd = (1/n)*(sum(singular_values)^2)

    ss_svd*P_eps_delta
end

"""
return the lower bound for the "Query andWorkload Error"
    according to the article:
    "Ecient Batch Query Answering Under Dierential Privacy"
"""
function getLowerBound(queries::RangeQueries, eps, delta)

    query_matrix = zeros(Int64, length(queries.intervals), queries.domain)

    for (index, interval) in enumerate(queries.intervals)
        query_matrix[index,:] = get_query_vector_from_interval(interval, queries.domain)
    end

    m,n = size(query_matrix)
    P_eps_delta = (2*log2(2/delta))/(eps^2)
    singular_values = svdfact(query_matrix)[:S]
    ss_svd = (1/n)*(sum(singular_values)^2)

    ss_svd*P_eps_delta
end
