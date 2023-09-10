using Plots;
using LinearAlgebra;

function thomas_algorithm(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64})
    n = length(d)
    c_prime = Vector{Float64}(undef, n - 1)
    d_prime = Vector{Float64}(undef, n)

    c_prime[1] = c[1] / b[1]
    d_prime[1] = d[1] / b[1]

    for i in 2:n-1
        c_prime[i] = c[i] / (b[i] - a[i-1] * c_prime[i-1])
    end

    for i in 2:n
        d_prime[i] = (d[i] - a[i-1] * d_prime[i-1]) / (b[i] - a[i-1] * c_prime[i-1])
    end

    x = Vector{Float64}(undef, n)
    x[n] = d_prime[n]

    for i in n-1:-1:1
        x[i] = d_prime[i] - c_prime[i] * x[i+1]
    end

    return x
end

function crank_nikolson(u₀::Function, λ::Float64, X::Vector{Float64}, T::Vector{Float64})
    u = zeros(length(X), length(T))

    for i in 1:length(X)
        u[i, 1] = u₀(X[i])
    end

    N = length(X)
    for n in 1:(length(T)-1)
        # bc
        u[1, n] = 0
        u[N, n] = 0

        A = zeros(N, N)
        B = zeros(N)

        A[1, 1] = 2 * (1 + λ)
        A[1, 2] = -λ
        B[1] = 2 * (1 - λ) * u[1, n] + λ * u[2, n]
        for i in 2:(N-1)
            A[i, i-1] = -λ
            A[i, i] = 2 * (1 + λ)
            A[i, i+1] = -λ
            B[i] = λ * u[i-1, n] + 2 * (1 - λ) * u[i, n] + λ * u[i+1, n]
        end
        A[N, N] = 2 * (1 + λ)
        A[N, N-1] = -λ
        B[N] = 2 * (1 - λ) * u[N, n] + λ * u[N-1, n]

        a = diag(A, -1)
        b = diag(A)
        c = diag(A, 1)
        # println(A)
        uᵢ = thomas_algorithm(a, b, c, vec(B))
        for i in 1:N
            u[i, n+1] = uᵢ[i]
        end
    end

    result = zeros(N)
    for i in 1:N
        result[i] = u[i, length(T)-1]
    end

    return result
end

function ParabolicNumerical(u₀, dx, F::{Float64}, timeStep)
    n = Int(1 / dx) + 1

    X = 0:dx:1

    u = zeros(n, timeStep)

    for i in 1:n
        u[i, 1] = u₀(X[i])
    end

    for i in 2:(timeStep)
        for j in 2:(n-1)
            u[j, i] += u[j, i-1] + (F / 2) * (u[j+1, i-1] - 2 * u[j, i-1] + u[j-1, i-1])
        end
        u[1, i] = u[1, i-1]
        u[n, i] = u[n, i-1]
    end

    res = zeros(n)
    for i in 1:n
        res[i] = u[i, timeStep-1]
    end

    return res
end


# (i)
α = 1 / (π^2)
Δt = 0.01
Δx = 0.1
X = 0:Δx:1
T = 0:Δt:0.1
u₀(x) = sin(π * x)
λ = (α * Δt) / (Δx^2)
uf = crank_nikolson(u₀, λ, Vector(X), Vector(T))
print(uf)
plot(X, uf, color=:red)

xspace = 0:(1/50):1
uexplicit = ParabolicNumerical(u₀, 1 / 50, 1, 11)
plot!(xspace, uexplicit, color=:blue)
png("crank_nikolson1.png")


# (ii)
α = 3
Δt = 0.01
Δx = 0.52
X = 0:Δx:π
T = 0:Δt:0.1
u₀(x) = sin(x)
λ = (α * Δt) / (Δx^2)
uf = crank_nikolson(u₀, λ, Vector(X), Vector(T))
print(uf)
plot(X, uf, color=:red)

xspace = 0:(π/50):π
uexplicit = ParabolicNumerical(u₀, π / 50, 1, 11)
plot!(xspace, uexplicit, color=:blue)
png("crank_nikolson2.png")