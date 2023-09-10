using Plots;
using LinearAlgebra;

function thomas_algorithm(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64})
    n = length(d)
    c_prime = Vector{Float64}(undef, n-1)
    d_prime = Vector{Float64}(undef, n)

    c_prime[1] = c[1] / b[1]
    d_prime[1] = d[1] / b[1]

    for i in 2:n-1
        c_prime[i] = c[i] / (b[i] - a[i-1]*c_prime[i-1])
    end

    for i in 2:n
        d_prime[i] = (d[i] - a[i-1]*d_prime[i-1]) / (b[i] - a[i-1]*c_prime[i-1])
    end

    x = Vector{Float64}(undef, n)
    x[n] = d_prime[n]

    for i in n-1:-1:1
        x[i] = d_prime[i] - c_prime[i]*x[i+1]
    end

    return x
end

function BVP2(h, l_x, r_x)
    n = (r_x - l_x)/h + 1
    n = Int(n)
    x = l_x:h:r_x
    result = zeros(n, 1)
    result[n] = 5
    coeff_matrix = zeros(n, n)
    coeff_matrix[1,1] = 1 
    coeff_matrix[n,n] = 1 

    for i in 2:n-1 
        coeff_matrix[i, i-1] = (1/(h^2)) - (1/(2*h))*(x[i]^2)
        coeff_matrix[i, i] = -2 * (2*x[i] + 1/(h^2))
        coeff_matrix[i, i+1] = (1 / (h^2)) + (1 / (2*h)) * (x[i]^2)
        result[i] = 0
    end

    # y_pred = (coeff_matrix^-1) * result
    a = diag(coeff_matrix, -1)
    b = diag(coeff_matrix)
    c = diag(coeff_matrix, 1)

    y = thomas_algorithm(a, b, c, vec(result))

    return y, x
end

y, x = BVP2(0.1, 0, 1)
println(y)
plot(x, y, color=:red)

sol(x) = x^4 + 4*x
plot!(0:0.1:1, sol.(0:0.1:1), color=:blue)

png("Ass5_2.png")