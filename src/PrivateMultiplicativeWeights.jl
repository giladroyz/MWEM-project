module PrivateMultiplicativeWeights

using
    Distributions: Laplace, wsample
using
    Printf,
    Hadamard,
    LinearAlgebra,
    Random,
    IterTools,#.subsets
    Statistics
    #Iterators.subsets

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
    kl_divergence_error,
    mean_squared_error,
    queriesMatrix,
    gosper,
    #get_parities,
    evaluate,
    get_parities,
#    getLowerBound,
    Interval,
    get_query_vector_from_interval

import
    Base: eltype, length, iterate
#Base: start, next, done, eltype, length

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
