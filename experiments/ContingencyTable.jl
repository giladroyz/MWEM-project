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

function run_full_test(contingencyTable::Histogram, epsilons::Array{Float64}, number_of_tests::Int,
    order_of_parities::Int, mwem_iterations::Int)

    @assert number_of_tests >= 1 "number of tests must be >= 1"
    @assert length(epsilons) >= 1 "must be at least one epsilon"
    @assert order_of_parities >= 1 "order of the parity queries must be >= 1"
    @assert mwem_iterations >= 1 "mwem iterations must be >= 1"

    number_of_samples = contingencyTable.num_samples

    queries = Parities(length(domain_dim), order_of_parities)

    mwem_errors = zeros(length(epsilons), number_of_tests) # the error vector of mwem
    #svd_errors = zeros(length(epsilons), number_of_tests) # the error vector of svd

    for test=1:number_of_tests

        println("test: ", test)

        for (e_index, e) in enumerate(epsilons)
            mwem_error = run_test(contingencyTable, e, queries, number_of_samples, mwem_iterations, false)
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


function main4()

    dimension = 8
    idx = [1,3,5,7]
    inx_out = []
    #bii = BinaryItr(dimension, idx)
    for s in BinaryItr(dimension, idx)
        push!(inx_out, s)
    end

    print(inx_out)

end

function norm_1_with_indices(alpha::Int64)

    @assert alpha >= 0

    count = 0
    index = 1
    idx = []
    while alpha > 0
        if alpha % 2 == 1
            count += 1
            push!(idx, index)
        end
        alpha >>= 1
        index += 1
    end
    count, Array{Int64}(idx)
end

function culc_marginal(queries::Parities, h::Histogram, alpha::Int64)

    alpha_norm1, alpha_indices = norm_1_with_indices(alpha)

    inverse_alpha = setdiff(collect(1:queries.dimension), alpha_indices)

    marginal = zeros(Float64, 2^alpha_norm1)

    for (beta_index, beta) in enumerate(BinaryItr(queries.dimension, alpha))
        beta_norm1, beta_indices = norm_1_with_indices(beta)
        temp_sum = 0
        for gamma in BinaryItr(queries.dimension, inverse_alpha, beta)
            temp_sum += h.weights[gamma+1]
        end
        marginal[beta_index] = temp_sum
    end

    marginal

end

function main5()

    nltcs = CSV.read("data/nltcs_new_csv_2.csv")
    dropmissing!(nltcs)

    #println(nltcs[1])

    contTable = readContingencyTable(nltcs[[:1, :2]], 16)

    number_of_samples = size(contTable)[1]

    hist = Histogram(contTable, number_of_samples)

    queries = Parities(16, 3)

    #print(hist.weights[1:10])
    #coeff = fourierCoefficients(queries, hist)
    
    time2 = @elapsed culc_marginal(queries, hist, 7)
    time1 = @elapsed calculateMarginal(queries, hist, 7)

    println(time1)
    println(time2)

end

main5()
