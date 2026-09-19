# ==============================================================================
# M466 Lab - System Information
#
# Machine:         Meridian1
# OS:              Windows 11 / Ubuntu WSL
#
# CPU:             Intel Core Ultra 5 245HX
# WSL CPUs:        14
# Architecture:    x86_64
# RAM:             31 GiB available to WSL
#
# GPU:             NVIDIA RTX PRO 1000 Blackwell Generation Laptop GPU
# Architecture:    NVIDIA Blackwell
# CUDA Cores:      2,560
# Tensor Cores:    80 (5th generation)
# RT Cores:        20 (4th generation)
# GPU Memory:      8 GB GDDR7 (~7.96 GiB usable)
# Memory Bus:      128-bit
# Memory Bandwidth: 384 GB/s
# Compute Cap.:    12.0 (sm_120)
# GPU Power Cap:   110 W
#
# NVIDIA Driver:   596.53
# Driver CUDA:     13.2
# CUDA.jl Runtime: 13.4.0
#
# Julia:           1.13.0
# LLVM:            20.1.8
# Julia Project:   ~/M466
#------------------------------------------------------------------------------
# Commands used to obtain system information
#
# WSL/Linux terminal:
#
#   lscpu
#       CPU model, architecture, and CPUs available to WSL
#
#   free -h
#       RAM and swap available to WSL
#
#   nvidia-smi
#       GPU model, VRAM, power limit, driver, and supported CUDA version
#
# Julia REPL:
#
#   versioninfo()
#       Julia, OS, CPU, LLVM, and thread information
#
#   using CUDA
#   dev = CUDA.device()
#   CUDA.capability(dev)
#       GPU device and compute capability
#
#   CUDA.versioninfo()
#       CUDA runtime, driver, libraries, Julia CUDA packages, and GPU information
#
# CUDA core count:
#   Obtained from published specifications for the exact GPU model;
#   nvidia-smi and CUDA.versioninfo() do not directly report CUDA core count.
# ==============================================================================



# flops.jl

using LinearAlgebra
using CUDA
using Plots

function main()

    warmpup = 100
    loop = [1000, 2000, 3000, 4000, 5000, 6000]
    trials = 5
    T = Float64
    cpu_Tmin = zeros(length(loop))
    
    println("n = ", warmpup, " warm-up iterations")
    println("precision = ", T)
    println()

    # --------------------------------------------------
    # Warm-up - Initialize GPU and compile kernels
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
    # CPU - LAPACK+BLAS
    # --------------------------------------------------

    println("CPU")
    BLAS.set_num_threads(7)  # note: CPU has 14 logical cores, but we limit BLAS to 7 threads.
                             # Setting BLAS threds to 7 had better performance than 14 threads on this CPU for this problem size.
    println("BLAS threads = ", BLAS.get_num_threads())

    cpu_flmax = 0.0

    for j in 1:length(loop)

        n = loop[j]
        println("Matrix size n = ", n)

        Tmin = Inf

        for i = 1:trials

            A = rand(T, n, n)
            xtrue = ones(T, n)
            b = A * xtrue

            Tn = @elapsed x = A \ b

            flop_count(n) = (2/3) * n^3 + 2 * n^2
            flops = flop_count(n) / Tn
            cpu_flmax = max(cpu_flmax, flops)

            error = norm(x - xtrue) / norm(xtrue)

            # Minimum time among the 5 trials
            Tmin = min(Tmin, Tn)

            println("Trial ", i)
            println("  Tn     = ", Tn)
            println("  Error  = ", error)
            println("  GFLOPS = ", flops / 1e9)
        end

    cpu_Tmin[j] = Tmin

end

p = plot(
    log.(loop),
    log.(cpu_Tmin),
    marker = :circle,
    legend = false,
    xlabel = "log(n)",
    ylabel = "log(Tn)",
    title = "CPU Execution Time vs Matrix Size"
)

savefig(p, "cpu_runtime.png")


end

main()