
module D1Q2Burgers

using LinearAlgebra: norm
using Printf: @printf

export d1q2, u0_ramp, burgers_exact, convergence

"Number of discrete velocities."
const Q = 2

"Velocity directions; the actual particle velocities are `λ .* C`."
const C = (-1, 1)


@inline feq(i::Int, u, λ) = u / 2 + C[i] * u^2 / (4λ)


function collision!(f::AbstractMatrix, s, λ)
    @inbounds for j in axes(f, 2)
        u = f[1, j] + f[2, j]
        for i in 1:Q
            f[i, j] += s * (feq(i, u, λ) - f[i, j])
        end
    end
    return f
end


function transport!(f::AbstractMatrix)
    Nx = size(f, 2)
    @inbounds for j in 1:(Nx - 1)     # c₁ = -1: f₁[j] ← f₁[j+1]
        f[1, j] = f[1, j + 1]
    end
    @inbounds for j in Nx:-1:2        # c₂ = +1: f₂[j] ← f₂[j-1]
        f[2, j] = f[2, j - 1]
    end
    return f
end


function d1q2(; k::Integer = 8, s::Real = 2.0, λ::Real = 1.0, T::Real = 0.2,
              u0 = u0_ramp)
    0 < s <= 2 || @warn "s = $s is outside the stability range (0, 2]"

    N  = 2^k
    dx = 1 / N
    dt = dx / λ
    nt = round(Int, T / dt)
    Nx = N + 1
    x  = range(0, 1; length = Nx)

    f = Matrix{Float64}(undef, Q, Nx)
    @inbounds for j in 1:Nx, i in 1:Q
        f[i, j] = feq(i, u0(x[j]), λ)
    end

    for _ in 1:nt
        collision!(f, s, λ)
        transport!(f)
    end

    u = [f[1, j] + f[2, j] for j in 1:Nx]
    return (; x = collect(x), u, k, s, λ, dx, dt, nt, T, t = nt * dt)
end

"""
    u0_ramp(x)

Initial datum: a C1, monotonically increasing ramp from -1 to +1, with a cubic
transition layer of half-width 1/4 centred at x = 1/2.

Being non-decreasing, it generates no shock: the characteristics never cross
and the entropy solution stays smooth for all time, which is what makes the
characteristics-based `burgers_exact` valid. Outside `[0, 1]` it is extended by
the constants ±1, consistently with the outflow boundary condition of the scheme.
"""
function u0_ramp(x)
    ξ = x - 0.5
    return abs(ξ) >= 0.25 ? sign(ξ) : sign(ξ) * (1 + (4 * abs(ξ) - 1)^3)
end


function burgers_exact(X, t; u0 = u0_ramp, lo = -1.0, hi = 2.0, iters::Integer = 80)
    a, b = float(lo), float(hi)
    for _ in 1:iters
        mid = 0.5 * (a + b)
        if mid + t * u0(mid) < X
            a = mid
        else
            b = mid
        end
    end
    return u0(0.5 * (a + b))
end


function convergence(; s_list = [2.0, 1.9, 1.75, 1.0, 0.75, 0.5],
                     ks = 3:16, λ::Real = 1.0, T::Real = 0.2, verbose::Bool = true)
    ks  = collect(ks)
    dxs = [2.0^(-k) for k in ks]

    errors = Dict{Float64, Vector{Float64}}()
    orders = Dict{Float64, Vector{Float64}}()

    for s in s_list
        err = Float64[]
        for k in ks
            sim  = d1q2(; k, s, λ, T)
            uex  = burgers_exact.(sim.x, sim.t)   # note: sim.t, not T
            push!(err, norm(sim.u .- uex) / norm(uex))
        end
        ord = [log(err[i] / err[i - 1]) / log(dxs[i] / dxs[i - 1]) for i in 2:length(ks)]
        errors[s], orders[s] = err, ord

        if verbose
            println("s = $s")
            for i in eachindex(ks)
                if i == 1
                    @printf("   k = %2d   dx = %.4e   err = %.4e\n", ks[i], dxs[i], err[i])
                else
                    @printf("   k = %2d   dx = %.4e   err = %.4e   order = %.2f\n",
                            ks[i], dxs[i], err[i], ord[i - 1])
                end
            end
        end
    end

    return (; ks, dxs, s_list, errors, orders)
end

end # module
