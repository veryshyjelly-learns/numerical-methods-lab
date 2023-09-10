using Plots;

function ParabolicNumerical(dx, F, timeStep)
    n = Int(1 / dx) + 1

    X = 0:dx:1

    u = zeros(n, timeStep)

    for i in 1:n
        u[i, 1] = sin(π * X[i])
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

xspace = 0:(1/50):1

u300 = ParabolicNumerical(1 / 50, 0.9, 300)
u400 = ParabolicNumerical(1 / 50, 0.9, 400)
u500 = ParabolicNumerical(1 / 50, 0.9, 500)

plot(xspace, u300, color=:red)
plot!(xspace, u400, color=:blue)
plot!(xspace, u500, color=:green)

png("parabolic1.png")