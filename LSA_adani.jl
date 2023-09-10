using Plots

function GaussElimination(A::Matrix{Float64}, b::Vector{Float64}, n::Int64)::Vector{Float64}
    for i in 1:n-1
        if A[i, i] == 0
            continue
        end
        for j in i+1:n
            m = A[j,i]/A[i,i]
            A[j,:] -= m * A[i,:]
            b[j] -= m * b[i]
        end
    end

    x = zeros(n)
    for i in n:-1:1
        x[i] = (b[i] - dot(A[i,i+1:n], x[i+1:n])) / A[i,i]
    end

    return vec(x)
end

function dot(A::Vector{Float64}, B::Vector{Float64})::Float64
    sum = 0.0
    for i in 1:length(A)
        sum += A[i] * B[i]
    end
    return sum
end

function LeastSquareAppr(x::Vector{Float64}, y::Vector{Float64},n::Int64)::Vector{Float64}
    m = length(x)
    A = zeros(n, n)
    b = zeros(n)

    for i in 1:n
        for j in 1:n
            for k in 1:m
                A[i,j] += x[k]^(2*n-i-j)
            end
        end
        for k in 1:m
            b[i] += y[k] * (x[k]^(n-i))
        end
    end

    # A/b

    return GaussElimination(A, b, n)
end

function Polynomial(coeff::Vector{Float64})::Function
    return x -> sum([coeff[i] * x^(length(coeff)-i) for i in 1:length(coeff)])
end

Ax::Vector{Float64} = vec([1, 1.5, 2, 2.5, 3, 3.5, 4])
Ay::Vector{Float64} = vec([25, 31, 27, 28, 36, 35, 32])
Aresult = LeastSquareAppr(Ax, Ay, 2)
println("1. y = ax + b, [a, b] = ", Aresult)
plot(1:0.5:4, Polynomial(Aresult).(1:0.5:4), label="y = ax + b", xlims=(0, 5), ylims=(20, 40))
scatter!(Ax, Ay, label="data")
png("1.png")

