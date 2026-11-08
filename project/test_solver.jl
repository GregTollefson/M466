
# test_solver.jl
# Phase 1: CPU multigrid convergence benchmark

using LinearAlgebra
using Printf

include("operator.jl")
include("redblack.jl")
include("residual.jl")
include("multigrid.jl")

function test_multigrid(N; tol=1e-8, maxcycles=50)

    h = 1.0 / (N + 1)

    u = zeros(Float64, N+2, N+2)
    f = ones(Float64, N, N)

    cycles = 0
    res = relative_residual(u, f, h)

    start_time = time_ns()

    while res > tol && cycles < maxcycles
        vcycle!(u, f, h)
        cycles += 1
        res = relative_residual(u, f, h)
    end

    elapsed = (time_ns() - start_time) / 1e9

    return cycles, res, elapsed
end

# Run convergence tests
grid_sizes = [3, 7, 15, 31, 63, 255]

println("\nCPU Multigrid Benchmark")
println("-"^65)

@printf("%8s %12s %12s %16s %12s\n",
        "N", "Unknowns", "V-cycles",
        "Rel. residual", "Time (s)")

println("-"^65)

for N in grid_sizes
    cycles, res, elapsed = test_multigrid(N)

    @printf("%8d %12d %12d %16.4e %12.6f\n",
            N, N^2, cycles, res, elapsed)
end

println("-"^65)
