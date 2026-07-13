# Convergence of the D1Q2 scheme as a function of the relaxation parameter s.
#   julia --project=examples examples/convergence.jl
# Expected: order 2 for s = 2 (the O(dx) term of the equivalent equation carries


using D1Q2Burgers
using Plots

s_list = [2.0, 1.9, 1.75, 1.0, 0.75, 0.5]
res = convergence(; s_list, ks = 3:16, T = 0.2)

plt = plot(; xscale = :log10, yscale = :log10,
           xlabel = "dx", ylabel = "relative L² error",
           title = "D1Q2 — Burgers", legend = :bottomright)

for s in s_list
    plot!(plt, res.dxs, res.errors[s]; lw = 2, marker = :circle, ms = 3, label = "s = $s")
end

# slope guides
plot!(plt, res.dxs, res.dxs;      ls = :dot, lc = :black, label = "O(dx)")
plot!(plt, res.dxs, res.dxs .^ 2; ls = :dash, lc = :black, label = "O(dx²)")

savefig(plt, joinpath(@__DIR__, "..", "convergence.png"))
display(plt)
