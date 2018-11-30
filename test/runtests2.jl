include("../src/PrivateMultiplicativeWeights.jl")

using PrivateMultiplicativeWeights
using Distributions
using Plots
#plotlyjs()
#plotly()
#pyplot()
gr()

d_continuous = Truncated(MixtureModel(Normal, [(.275, .075), (.65, .1)], [0.65, 0.35]), 0, 1)
xs = linspace(d_continuous.lower, d_continuous.upper, 300)
ys = pdf.(d_continuous, xs)
plot(xs, ys);

num_samples = 1000
domain_size = 100
samples = rand(d_continuous, num_samples)

function histogram_from_samples(samples, num_bins)
    @assert 0 <= minimum(samples)
    @assert maximum(samples) <= 1
    h = zeros(num_bins)
    new_samples = samples*num_bins
    for x in round.(Int64, new_samples)
        h[x] += 1
    end
    collect(1:num_bins)/num_bins, h
end

xs, ys = histogram_from_samples(samples, domain_size)
bar(xs, ys, width=1/domain_size);
#gui()

mw = mwem(RangeQueries(domain_size), Histogram(ys/sum(ys), num_samples))
bar(xs, mw.synthetic.weights*num_samples, color="red", alpha=0.3, width=1/domain_size);
#gui()

mw = mwem(RangeQueries(domain_size), Histogram(ys/sum(ys), num_samples), MWParameters(epsilon=100, iterations=80))
bar(xs, mw.synthetic.weights*num_samples, color="red", alpha=0.3, width=1/domain_size);
#gui()

a = bar(xs, ys)#, alpha=1/domain_size))
bar!(xs, mw.synthetic.weights*num_samples, color="red", alpha=0.3, width=1/domain_size)

display(a)
