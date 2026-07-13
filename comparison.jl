# Numerical solution vs. exact solution of the Burgers equation.
#
#   julia --project=examples examples/comparison.jl

using D1Q2Burgers
using Plots

res  = d1q2(; k = 10, s = 2.0, T = 0.2)
u_ex = burgers_exact.(res.x, res.t)     # compare at the time actually reached

plt = plot(res.x, res.u; lw = 2, label = "D1Q2 (s = $(res.s))",
           xlabel = "x", ylabel = "u", title = "Burgers, t = $(round(res.t, digits = 4))")
plot!(plt, res.x, u_ex; lw = 1, ls = :dash, label = "exact")

savefig(plt, joinpath(@__DIR__, "..", "comparison.png"))
display(plt)
