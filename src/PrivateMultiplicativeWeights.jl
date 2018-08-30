module PrivateMultiplicativeWeights

using
    Distributions.Laplace,
    Distributions.wsample,
    Hadamard,
    Iterators.subsets

export
    mwem,
    MWParameters,
    Tabular,
    Histogram,
    HistogramQueries,
    SeriesRangeQueries,
    RangeQueries,
    RangeQuery,
    Parities,
    FactorParities,
    maximum_error,
    mean_squared_error,
    queriesMatrix,
#    getLowerBound,
    Interval,
    get_query_vector_from_interval

import
    Base: start, next, done, eltype, length

include("interface.jl")
include("histogram.jl")
include("rangequeries.jl")
include("factors.jl")
include("gosper.jl")
include("parities.jl")
include("error.jl")
include("mwem.jl")
#include("MatrixMechanism.jl")

"""
    PrivateMultiplicativeWeights

A simple and practical algorithm for differentially private data analysis.
"""
PrivateMultiplicativeWeights

end # module
