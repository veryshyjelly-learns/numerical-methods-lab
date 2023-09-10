
mutable struct Point
    y::Float64
    x1::Float64
    x2::Float64
end

function NewtonInterpolation(x::Vector{Float64}, y::Vector{Float64}, n::Int)::Vector{Float64}
    differences = Vector{Vector{Point}}()

    push!(differences, Vector{Point}())
    for i in 1:n
        push!(differences[1], Point(y[i], x[i], x[i]))
    end

    idx = 1
    for _ in 2:n
        push!(differences, Vector{Point}())
        for i in 2:length(differences[idx])
            firstpoint::Point = differences[idx][i-1]
            secondpiont::Point = differences[idx][i]

            X2::Float64 = secondpiont.x2 
            X1::Float64 = firstpoint.x1 
            Y::Float64 = (secondpiont.y - firstpoint.y) / (X2 - X1)

            push!(differences[idx+1], Point(Y, X1, X2))
        end
    end

    coeffs = Vector{Float64}()
    for diff in differences 
        push!(coeffs, diff[1].y)
    end

    return coeffs
end

print("Enter the number of points: ")
n = parse(Int, readline())

x = Vector{Float64}()
y = Vector{Float64}()

println("Enter the values of x")
valuesX = map(a -> parse(Float64, a), split(strip(readline()), " "))

for i in 1:n 
    push!(x, valuesX[i])
end


println("Enter the values of y")
valuesY = map(a -> parse(Float64, a), split(strip(readline()), " "))

for i in 1:n 
    push!(y, valuesY[i])
end

coefficients = NewtonInterpolation(x, y, n)
println(coefficients)
