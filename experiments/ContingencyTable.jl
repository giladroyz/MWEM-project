include("DataAccessibility.jl")
push!(LOAD_PATH, abspath("./src/"))
#include(abspath("./src/parities.jl"))

using PrivateMultiplicativeWeights
using Distributions
using Plots
using DataFrames
using CSV

gr()

function run_test(data::Histogram, epsilon::Float64, queries::Parities,
    number_of_samples::Int, number_iterations::Int, noisy_init::Bool)

    mw = mwem(queries, data,
        MWParameters(epsilon=epsilon, iterations=number_iterations, noisy_init=noisy_init))
        #println(sum(mw.real.weights))
        #println(sum(mw.synthetic.weights))
    mwem_error = kl_divergence_error(mw)#mean_squared_error(mw)*(number_of_samples^2)
    
    mwem_error
end

function run_full_test(data::DataFrames.DataFrame, epsilons::Array{Float64}, number_of_tests::Int,
    order_of_parities::Int, mwem_iterations::Int)

    @assert number_of_tests >= 1 "number of tests must be >= 1"
    @assert length(epsilons) >= 1 "must be at least one epsilon"
    @assert order_of_parities >= 1 "order of the parity queries must be >= 1"
    @assert mwem_iterations >= 1 "mwem iterations must be >= 1"

    number_of_samples = size(data)[1]

    # Tuple with length = number of columns, and each entry has the number of elements in the column
    domain_dim = getDomainDim2(data)

    # the normalized data in a matrix
    data_matrix = covertData(data)

    #println(size(data_matrix))

    barak_alg_data = Tabular(data_matrix)
    data_histogram = Histogram(barak_alg_data)

    println(length(domain_dim))

    queries = Parities(length(domain_dim), order_of_parities)

    mwem_errors = zeros(length(epsilons), number_of_tests) # the error vector of mwem
    #svd_errors = zeros(length(epsilons), number_of_tests) # the error vector of svd

    for test=1:number_of_tests

        println("test: ", test)

        for (e_index, e) in enumerate(epsilons)
            mwem_error = run_test(data_histogram, e, queries, number_of_samples, mwem_iterations, false)
            mwem_errors[e_index, test] = mwem_error
        end
    end

    mwem_errors#, svd_errors

end

function write_results(file_name::String, mwem_errors::AbstractArray, epsilons::AbstractArray)

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

    close(io)

end

function main()

    ### load datasets ###
    nltc = CSV.read("data/nltcs.csv")
    #####################

    ### drop rows with missing values
    dropmissing!(nltc)
    #####################

    epsilons = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.06, 0.07, 0.08, 0.09, 0.1]
    order = 3
    mwem_iterations = 10
    number_of_tests = 100

    ## nltc
    println("nltc")
    dataset = nltc
    mwem_errors = run_full_test(dataset, epsilons, number_of_tests, order, mwem_iterations)
    write_results("nltc_result.txt", mwem_errors, epsilons)

end

function main2()
    p = Parities(4,1)
    #g = [s for s in gosper(10, 3)]
    println(p.idx)
    println(length(p.idx))
    for i in 1:5
        #a = get(p, i)
        #println(a.weights)
        #println(length(a.weights))
    end
    #println(g)
    #println(length(g))
end

function main3()

    nltc = CSV.read("data/nltcs.csv")
    #####################

    ### drop rows with missing values
    dropmissing!(nltc)
    
    data = nltc[[:1, :2, :3, :4, :5, :6]]

    epsilons = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.06, 0.07, 0.08, 0.09, 0.1]
    order = 3
    mwem_iterations = 10
    number_of_tests = 10

    domain_dim = getDomainDim2(data)

    data_matrix = covertData(data)

    barak_alg_data = Tabular(data_matrix)
    data_histogram = Histogram(barak_alg_data)

    queries = Parities(length(domain_dim), order)

    mw = mwem(queries, data_histogram, MWParameters(epsilon=epsilons[1], iterations=mwem_iterations, noisy_init=false))

    # io = open("nltcs_histogram.txt", "w")

    # for a in  mw.synthetic.weights*size(data_matrix)[1]
    #     println(io, string(a))
    # end
    # close(io)

    # io = open("nltcs_histogram_real.txt", "w")

    # for a in  mw.real.weights*size(data_matrix)[1]
    #     println(io, string(a))
    # end
    # close(io)

    println(evaluate(queries, mw.real))
    println(evaluate(queries, mw.synthetic))
    println(length(evaluate(queries, mw.synthetic)))
    println(get_parities(queries,2))
    println(evaluate(get_parities(queries,2), mw.synthetic))
    #bar([1:2^16], mw.synthetic.weights*size(data_matrix)[1], color="red", alpha=0.3, width=1/2^16)

end

main()
