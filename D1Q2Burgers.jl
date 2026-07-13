"""
    D1Q2Burgers

A minimal D1Q2 lattice Boltzmann scheme for the inviscid Burgers equation

    ∂ₜu + ∂ₓ(u²/2) = 0,   x ∈ [0, 1],

written as a single-relaxation-time (BGK) scheme, following Graille (2014).

See `d1q2` for the solver, `burgers_exact` for the reference solution and
`convergence` for the mesh-refinement study.
"""
module D1Q2Burgers

using LinearAlgebra: norm
using Printf: @printf

export d1q2, u0_ramp, burgers_exact, convergence

"Number of discrete velocities."
const Q = 2

"Velocity directions; the actual particle velocities are `λ .* C`."
const C = (-1, 1)

"""
    feq(i, u, λ)

Equilibrium distribution of the `i`-th population.

The two moments are `m₀ = f₁ + f₂ = u` (conserved) and `m₁ = λ(f₂ - f₁)`, whose
equilibrium is the Burgers flux `m₁ᵉᵠ = u²/2`. Inverting the moment matrix gives

    fᵢᵉᵠ(u) = u/2 + cᵢ u² / (4λ).
"""
@inline feq(i::Int, u, λ) = u / 2 + C[i] * u^2 / (4λ)

"""
    collision!(f, s, λ)

BGK relaxation, in place: `fᵢ ← fᵢ + s (fᵢᵉᵠ - fᵢ)`, where `s = 1/τ ∈ (0, 2]`.

Since `Σᵢ fᵢᵉᵠ(u) = u`, the density `m₀` is conserved exactly; only `m₁` relaxes.
"""
function collision!(f::AbstractMatrix, s, λ)
    @inbounds for j in axes(f, 2)
        u = f[1, j] + f[2, j]
        for i in 1:Q
            f[i, j] += s * (feq(i, u, λ) - f[i, j])
        end
    end
    return f
end

"""
    transport!(f)

Free flight, in place: each population is shifted by exactly one cell
(`λ dt = dx`).

The incoming populations at the two ends (`f₁` at the right boundary, `f₂` at the
left one) are left untouched, which amounts to a zero-gradient / outflow boundary
condition. This is harmless for the test case of `u0_ramp`, whose solution stays
constant (±1) near the boundaries over the simulated time.
"""
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

"""
    d1q2(; k=8, s=2.0, λ=1.0, T=0.2, u0=u0_ramp)

Run the D1Q2 scheme on `[0, 1]` with `N = 2^k` cells, relaxation parameter `s`,
lattice velocity `λ` and final time `T`.

The acoustic scaling `dt = dx / λ` is used, so the number of time steps is
`nt = round(Int, T/dt)` and the time actually reached, `t = nt*dt`, generally
differs slightly from `T`. **Always compare against the exact solution at `t`,
not at `T`.**

Stability requires `0 < s ≤ 2` and the sub-characteristic condition
`λ ≥ max|u|` (here `max|u| = 1`, so `λ = 1` is the marginal choice).

Returns a `NamedTuple` `(; x, u, k, s, λ, dx, dt, nt, T, t)`.
"""
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

Initial datum: a C¹, monotonically increasing ramp from -1 to +1, with a cubic
transition layer of half-width 1/4 centred at `x = 1/2`.

Being non-decreasing, it generates **no shock**: the characteristics never cross
and the entropy solution stays smooth for all time, which is what makes the
characteristics-based `burgers_exact` valid. Outside `[0, 1]` it is extended by
the constants ±1, consistently with the outflow boundary condition of the scheme.
"""
function u0_ramp(x)
    ξ = x - 0.5
    return abs(ξ) >= 0.25 ? sign(ξ) : sign(ξ) * (1 + (4 * abs(ξ) - 1)^3)
end

"""
    burgers_exact(X, t; u0=u0_ramp, lo=-1.0, hi=2.0, iters=80)

Exact solution of the inviscid Burgers equation at `(X, t)`, obtained by tracing
the characteristic back to its foot `x₀`, i.e. by solving

    x₀ + t u₀(x₀) = X

with a bisection on `[lo, hi]`. The map `x₀ ↦ x₀ + t u₀(x₀)` is increasing (as
`u₀` is), so the root is unique and the bisection is guaranteed to converge.

`u(X, t) = u₀(x₀)`.
"""
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

"""
    convergence(; s_list, ks, λ=1.0, T=0.2, verbose=true)

Mesh-refinement study: for every relaxation parameter in `s_list`, compute the
relative discrete L² error against the exact solution for every `k` in `ks`, and
estimate the observed order between two consecutive meshes.

Returns `(; ks, dxs, s_list, errors, orders)` where `errors` and `orders` are
`Dict`s indexed by `s`.
"""
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
