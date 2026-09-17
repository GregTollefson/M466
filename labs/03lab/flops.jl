# flops.jl

using LinearAlgebra
using CUDA

function main()

    n = 20000
    
    trials = 5
    T = Float64

    flop_count(n) = (2/3) * n^3 + 2 * n^2

    println("n = ", n)
    println("precision = ", T)
    println()

    # --------------------------------------------------
    # Warm-up
    # --------------------------------------------------

    nw = 100

    A = rand(T, nw, nw)
    b = A * ones(T, nw)
    A \ b

    Ag = CuArray(A)
    bg = CuArray(b)

    CUDA.@sync begin
        Ag \ bg
    end

    # --------------------------------------------------
    # CPU
    # --------------------------------------------------

    println("CPU")
    println("BLAS threads = ", BLAS.get_num_threads())

    cpu_flmax = 0.0

    for i = 1:trials

        A = rand(T, n, n)
        xtrue = ones(T, n)
        b = A * xtrue

        Tn = @elapsed x = A \ b

        flops = flop_count(n) / Tn
        cpu_flmax = max(cpu_flmax, flops)

        error = norm(x - xtrue) / norm(xtrue)

        println("Trial ", i)
        println("  Tn     = ", Tn)
        println("  Error  = ", error)
        println("  GFLOPS = ", flops / 1e9)
    end

    println("Maximum CPU GFLOPS = ", cpu_flmax / 1e9)

    # --------------------------------------------------
    # GPU
    # --------------------------------------------------

    println()
    println("GPU")

    gpu_flmax = 0.0

    for i = 1:trials

        A = rand(T, n, n)
        xtrue = ones(T, n)
        b = A * xtrue

        Ag = CuArray(A)
        bg = CuArray(b)
        xtrue_g = CuArray(xtrue)

        CUDA.synchronize()

        Tn = @elapsed begin
            xg = Ag \ bg
            CUDA.synchronize()
        end

        flops = flop_count(n) / Tn
        gpu_flmax = max(gpu_flmax, flops)

        error = norm(xg - xtrue_g) / norm(xtrue_g)

        println("Trial ", i)
        println("  Tn     = ", Tn)
        println("  Error  = ", error)
        println("  GFLOPS = ", flops / 1e9)
    end

    println("Maximum GPU GFLOPS = ", gpu_flmax / 1e9)
end

main()