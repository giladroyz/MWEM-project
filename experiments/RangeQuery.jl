include("DataAccessibility.jl")
include("../src/PrivateMultiplicativeWeights.jl")

using PrivateMultiplicativeWeights
using Distributions
using Plots
using DataFrames
using CSV

include("../src/MatrixMechanism.jl")

#plotlyjs()
#plotly()
gr()

###################

"""
generate list of random range queries.
Params:
    number_of_queries: the number of range queries
    domain_size: the size of the domain.
    (if the domain is {1,...,N} than a range query is an interval I=(i,j) where 1<=i<=j<=N)
"""
function generate_random_range_queries(number_of_queries::Int, domain_size::Int)

    vector_range_queries = Array{Interval}(number_of_queries)
    for i=1:number_of_queries
        d = Categorical([1/domain_size for j=1:domain_size])
        a = rand(d,2)
        interval = (minimum(a), maximum(a))
        vector_range_queries[i] = interval
    end
    range_queries = RangeQueries(vector_range_queries, domain_size)

end

"""
run a test that compare the error of the mwem algorithm with the
svd lower bound of the matrix mechanism for the given range queries.

Params:
    ys: the y axis of the histogram that represent the data set.
    epsilon: the epsilon value of the e-privacy
    number_of_queries: the number of range queries that we want to give to the mwem algorithm
    number_of_samples: the number of rows in the dataset
    number_iterations: the number of iterations the mwem will do.
    noisy_init: how to init the mwem algorithm (noisy uniform, or uniform)
"""
function run_test(ys::AbstractArray, epsilon::Float64, range_queries::RangeQueries,
    number_of_samples::Int, number_iterations::Int, noisy_init::Bool)

    mw = mwem(range_queries, Histogram(ys/sum(ys), number_of_samples),
        MWParameters(epsilon=e, iterations=number_iterations, noisy_init=noisy_init))
    mwem_error = mean_squared_error(mw)
    svd_error = getLowerBound(range_queries, epsilon, 1/number_of_samples)

    (mwem_error, svd_error)
end

"""
run a full test that compare the error of the mwem algorithm with the
svd lower bound of the matrix mechanism for the given range queries.
it runs number of independent tests and give all the errors.

Params:
    data: the dataset fot the mwem to generate synthetic data
    epsilon: the epsilon value of the e-privacy
    number_of_queries: the number of range queries that we want to give to the mwem algorithm
    number_of_samples: the number of rows in the dataset
    number_iterations: the number of iterations the mwem will do.
    noisy_init: how to init the mwem algorithm (noisy uniform, or uniform)
"""
function run_full_test(data::DataFrames.DataFrame, epsilons::Array{Float64}, number_of_tests::Int,
    number_of_queries::Int, mwem_iterations::Int)

    @assert number_of_tests >= 1 "number of tests must be >= 1"
    @assert length(epsilons) >= 1 "must be at least one epsilon"
    @assert number_of_queries >= 1 "number of queries must be >= 1"
    @assert mwem_iterations >= 1 "mwem iterations must be >= 1"

    number_of_samples = size(data)[1]

    # Tuple with length = number of columns, and eaxh entry has the number of elements in the column
    domain_dim = getDomainDim2(data)

    # the normalized data in a matrix
    data_matrix = covertData(data)

    # the 1-D histogram of the elements of the domain in the dataset
    xs, ys = histogram_from_sample(data_matrix, domain_dim)

    mwem_errors = zeros(length(epsilons), number_of_tests) # the error vector of mwem
    svd_errors = zeros(length(epsilons), number_of_tests) # the error vector of svd

    for test=1:number_of_tests

        println("test: ", test)

        range_queries = generate_random_range_queries(number_of_queries, prod(domain_dim))

        for (e_index, e) in enumerate(epsilons)
            mwem_error, svd_error = run_test(ys, e, range_queries, number_of_samples, mwem_iterations, false)
            mwem_errors[e_index, test] = mwem_error
            svd_errors[e_index, test] = svd_error
        end
    end

    mwem_errors, svd_errors

end

function write_results(file_name::String, mwem_errors::AbstractArray,
    svd_errors::AbstractArray, epsilons::AbstractArray)

    num_epsilons, num_tests = size(mwem_errors)

    @assert length(epsilons) == num_epsilons "wrong number of epsilons"

    io = open(file_name, "w")

    println(io, num_epsilons)
    println(io, num_tests)
    print(io, "\n")

    for i=1:num_epsilons
        print(io, epsilons[i], ", ")
        epsilons
    end

    print(io, "\n")

    for i=1:num_tests
        for a in mwem_errors[:,i]
            print(io, a, ", ")
        end
        print(io, "\n")
    end

    print(io, "\n")

    for i=1:num_tests
        for a in svd_errors[:,i]
            print(io, a, ", ")
        end
        print(io, "\n")
    end

    close(io)

end

###################

function main()

    ### load datasets ###
    adult = CSV.read("data/adult.csv")
    transfusion = CSV.read("data/transfusion.csv")
    #####################

    ### drop rows with missing values
    dropmissing!(adult)
    dropmissing!(transfusion)
    #####################

    epsilons = [0.0125, 0.025, 0.1, 0.5]
    number_of_queries = 2000
    mwem_iterations = 10
    number_of_tests = 5

    ## Adult: age x hourss
    println("Adult: age x hourss")
    dataset = adult[[:1, :13]]
    (mwem_errors, svd_errors) = run_full_test(dataset, epsilons, number_of_tests, number_of_queries, mwem_iterations)
    write_results("Adult_age_hour.txt", mwem_errors, svd_errors, epsilons)

    ## Adult: capital loss
    # println("Adult: capital loss")
    # dataset =adult[[:12]]
    # println(getDomainDim2(dataset))
    # (mwem_errors, svd_errors) = run_full_test(dataset, epsilons, number_of_tests, number_of_queries, mwem_iterations)
    # write_results("Adult_capitalloss.txt", mwem_errors, svd_errors, epsilons)

    ## blood: recency x frequency
    println("blood: recency x frequency")
    dataset =transfusion[[:1, :2]]
    (mwem_errors, svd_errors) = run_full_test(dataset, epsilons, number_of_tests, number_of_queries, mwem_iterations)
    write_results("Transfusion_recency_frequency.txt", mwem_errors, svd_errors, epsilons)

    ### blood: monetary
    println("blood: monetary")
    dataset =transfusion[[:3]]
    (mwem_errors, svd_errors) = run_full_test(dataset, epsilons, number_of_tests, number_of_queries, mwem_iterations)
    write_results("Transfusion_monetary.txt", mwem_errors, svd_errors, epsilons)

end

main()
