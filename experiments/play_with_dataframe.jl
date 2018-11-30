using DataFrames

df = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
n, m = size(df)
println(n,m)
#println(df)
#println(df[1,:])

for i =1:m
    println(df[1,:][i][1])
end
#df2 = DataFrame([1,"L"])
push!(df, [1,"F"])
println(df)

a = DataFrame()
for i=1:5
    a[i] = []
end
push!(a, [1, "F", 2, 3, 5])
print(a)